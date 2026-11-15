"""SQLite 存储：报告与位点结果。

每次比对（无论成功或失败）在单个事务内写入自己的报告行与位点行，
不同批次的数据通过 report_id 严格隔离；旧报告只读，永不改写。
"""
from __future__ import annotations

import hashlib
import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

from .compare import ComparedSite

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
