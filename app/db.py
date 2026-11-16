"""SQLite 存储：报告与位点结果、坐标迁移结果。

每次比对（无论成功或失败）在单个事务内写入自己的报告行与位点行，
不同批次的数据通过 report_id 严格隔离；旧报告只读，永不改写。

迁移结果同样单事务写入（migrations + migrated_variants）；
migrations.dedup_key 唯一约束保证同一源报告、同一目标序列、同一映射内容
的重复提交（包括并发提交）只产生一份迁移报告。
"""
from __future__ import annotations

import hashlib
import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

from .compare import ComparedSite
from .migrate import MigratedRow

SCHEMA = """
CREATE TABLE IF NOT EXISTS reports (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    created_at TEXT NOT NULL,
    status TEXT NOT NULL,              -- 'done' | 'failed'
    error TEXT,                        -- 失败原因（含文件名与行号）
    ref_name TEXT,
    ref_length INTEGER,
    ref_sha256 TEXT,
    fasta_filename TEXT,
    vcf_a_filename TEXT,
    vcf_b_filename TEXT
);
CREATE TABLE IF NOT EXISTS variants (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    report_id INTEGER NOT NULL REFERENCES reports(id),
    pos INTEGER NOT NULL,
    ref TEXT NOT NULL,                 -- 归一化 REF，插入为 ''
    alt TEXT NOT NULL,                 -- 归一化 ALT，缺失为 ''
    category TEXT NOT NULL,            -- 'common' | 'a_only' | 'b_only'
    lines_a TEXT NOT NULL,             -- JSON 数组，甲 VCF 原行号
    lines_b TEXT NOT NULL              -- JSON 数组，乙 VCF 原行号
);
CREATE INDEX IF NOT EXISTS idx_variants_report ON variants(report_id);
CREATE TABLE IF NOT EXISTS migrations (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    created_at TEXT NOT NULL,
    source_report_id INTEGER NOT NULL REFERENCES reports(id),
    old_ref_name TEXT NOT NULL,
    old_ref_length INTEGER NOT NULL,
    old_ref_sha256 TEXT NOT NULL,
    new_ref_name TEXT NOT NULL,
    new_ref_length INTEGER NOT NULL,
    new_ref_sha256 TEXT NOT NULL,
    n_blocks INTEGER NOT NULL,         -- 映射块数
    mapping_sha256 TEXT NOT NULL,      -- 映射内容摘要（规范化块文本）
    dedup_key TEXT NOT NULL UNIQUE     -- 源报告:目标序列摘要:映射摘要
);
CREATE TABLE IF NOT EXISTS migrated_variants (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    migration_id INTEGER NOT NULL REFERENCES migrations(id),
    src_pos INTEGER NOT NULL,          -- 旧参考坐标（1-based）
    pos INTEGER NOT NULL,              -- 新参考坐标（1-based）
    ref TEXT NOT NULL,                 -- 迁移后 REF（负向块已取反向互补）
    alt TEXT NOT NULL,                 -- 迁移后 ALT
    category TEXT NOT NULL,            -- 沿用源报告归类
    lines_a TEXT NOT NULL,             -- JSON 数组，甲 VCF 原行号
    lines_b TEXT NOT NULL              -- JSON 数组，乙 VCF 原行号
);
CREATE INDEX IF NOT EXISTS idx_migrated_variants_migration
    ON migrated_variants(migration_id);
"""


def connect(db_path: str) -> sqlite3.Connection:
    # check_same_thread=False：FastAPI 依赖的建立/清理可能跑在不同工作线程，
    # 每个请求各自持有独立连接，请求内串行使用，是安全的。
    conn = sqlite3.connect(db_path, timeout=10, check_same_thread=False)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA busy_timeout = 10000")
    return conn


def init_db(db_path: str) -> None:
    Path(db_path).parent.mkdir(parents=True, exist_ok=True)
    conn = connect(db_path)
    try:
        conn.execute("PRAGMA journal_mode = WAL")  # 读写并发更稳
        conn.executescript(SCHEMA)
        conn.commit()
    finally:
        conn.close()


def ref_digest(ref_seq: str) -> str:
    return hashlib.sha256(ref_seq.encode("ascii")).hexdigest()


def insert_report(
    conn: sqlite3.Connection,
    *,
    status: str,
    error: str | None,
    ref_name: str | None,
    ref_length: int | None,
    ref_sha256: str | None,
    fasta_filename: str | None,
    vcf_a_filename: str | None,
    vcf_b_filename: str | None,
    compared: list[ComparedSite],
) -> int:
    """单事务写入一份报告及其全部位点；失败整体回滚，不留半成品。"""
    with conn:  # 成功提交、异常回滚
        cur = conn.execute(
            """INSERT INTO reports
               (created_at, status, error, ref_name, ref_length, ref_sha256,
                fasta_filename, vcf_a_filename, vcf_b_filename)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            (
                datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
                status,
                error,
                ref_name,
                ref_length,
                ref_sha256,
                fasta_filename,
                vcf_a_filename,
                vcf_b_filename,
            ),
        )
        report_id = cur.lastrowid
        conn.executemany(
            """INSERT INTO variants
               (report_id, pos, ref, alt, category, lines_a, lines_b)
               VALUES (?, ?, ?, ?, ?, ?, ?)""",
            [
                (
                    report_id,
                    row.site.pos,
                    row.site.ref,
                    row.site.alt,
                    row.category,
                    json.dumps(list(row.lines_a)),
                    json.dumps(list(row.lines_b)),
                )
                for row in compared
            ],
        )
    return report_id


def list_reports(conn: sqlite3.Connection) -> list[sqlite3.Row]:
    return conn.execute(
        """SELECT r.*,
                  (SELECT COUNT(*) FROM variants v
                    WHERE v.report_id = r.id AND v.category = 'common') AS n_common,
                  (SELECT COUNT(*) FROM variants v
                    WHERE v.report_id = r.id AND v.category = 'a_only') AS n_a_only,
                  (SELECT COUNT(*) FROM variants v
                    WHERE v.report_id = r.id AND v.category = 'b_only') AS n_b_only
           FROM reports r ORDER BY r.id DESC"""
    ).fetchall()


def get_report(
    conn: sqlite3.Connection, report_id: int
) -> tuple[sqlite3.Row | None, list[sqlite3.Row]]:
    report = conn.execute(
        "SELECT * FROM reports WHERE id = ?", (report_id,)
    ).fetchone()
    if report is None:
        return None, []
    variants = conn.execute(
        "SELECT * FROM variants WHERE report_id = ? ORDER BY pos, ref, alt",
        (report_id,),
    ).fetchall()
    return report, variants


# ---------- 坐标迁移 ----------


def list_migratable_reports(conn: sqlite3.Connection) -> list[sqlite3.Row]:
    """已完成、含 SNV 且无插入缺失的报告（迁移页可选的源报告）。"""
    return conn.execute(
        """SELECT r.*,
                  (SELECT COUNT(*) FROM variants v
                    WHERE v.report_id = r.id
                      AND length(v.ref) = 1 AND length(v.alt) = 1) AS n_snv
           FROM reports r
           WHERE r.status = 'done'
             AND EXISTS (SELECT 1 FROM variants v
                          WHERE v.report_id = r.id
                            AND length(v.ref) = 1 AND length(v.alt) = 1)
             AND NOT EXISTS (SELECT 1 FROM variants v
                              WHERE v.report_id = r.id
                                AND NOT (length(v.ref) = 1 AND length(v.alt) = 1))
           ORDER BY r.id DESC"""
    ).fetchall()


def insert_migration(
    conn: sqlite3.Connection,
    *,
    source_report_id: int,
    old_ref_name: str,
    old_ref_length: int,
    old_ref_sha256: str,
    new_ref_name: str,
    new_ref_length: int,
    new_ref_sha256: str,
    n_blocks: int,
    mapping_sha256: str,
    rows: list[MigratedRow],
) -> tuple[int, bool]:
    """单事务写入迁移报告及全部迁移位点，返回 (migration_id, 是否新建)。

    dedup_key（源报告 + 目标序列摘要 + 映射摘要）有唯一约束：
    内容相同的重复提交（含并发）由 INSERT OR IGNORE 去重，返回已有报告 id，
    不会产生重复或半份迁移报告。
    """
    dedup_key = f"{source_report_id}:{new_ref_sha256}:{mapping_sha256}"
    with conn:  # 成功提交、异常回滚
        cur = conn.execute(
            """INSERT OR IGNORE INTO migrations
               (created_at, source_report_id,
                old_ref_name, old_ref_length, old_ref_sha256,
                new_ref_name, new_ref_length, new_ref_sha256,
                n_blocks, mapping_sha256, dedup_key)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            (
                datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
                source_report_id,
                old_ref_name,
                old_ref_length,
                old_ref_sha256,
                new_ref_name,
                new_ref_length,
                new_ref_sha256,
                n_blocks,
                mapping_sha256,
                dedup_key,
            ),
        )
        if cur.rowcount == 1:
            migration_id = cur.lastrowid
            conn.executemany(
                """INSERT INTO migrated_variants
                   (migration_id, src_pos, pos, ref, alt, category, lines_a, lines_b)
                   VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
                [
                    (
                        migration_id,
                        row.site.src_pos,
                        row.site.pos,
                        row.site.ref,
                        row.site.alt,
                        row.category,
                        json.dumps(list(row.lines_a)),
                        json.dumps(list(row.lines_b)),
                    )
                    for row in rows
                ],
            )
            return migration_id, True
        # 内容相同的迁移已存在（可能由并发请求刚提交）：返回同一份结果
        row = conn.execute(
            "SELECT id FROM migrations WHERE dedup_key = ?", (dedup_key,)
        ).fetchone()
        return row["id"], False


def list_migrations(conn: sqlite3.Connection) -> list[sqlite3.Row]:
    return conn.execute(
        """SELECT m.*,
                  (SELECT COUNT(*) FROM migrated_variants v
                    WHERE v.migration_id = m.id) AS n_sites
           FROM migrations m ORDER BY m.id DESC"""
    ).fetchall()


def get_migration(
    conn: sqlite3.Connection, migration_id: int
) -> tuple[sqlite3.Row | None, list[sqlite3.Row]]:
    migration = conn.execute(
        "SELECT * FROM migrations WHERE id = ?", (migration_id,)
    ).fetchone()
    if migration is None:
        return None, []
    variants = conn.execute(
        "SELECT * FROM migrated_variants WHERE migration_id = ? ORDER BY pos, ref, alt",
        (migration_id,),
    ).fetchall()
    return migration, variants
