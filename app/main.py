"""FastAPI 网页服务：上传参考序列与甲乙两份 VCF，查看位点差异报告；
坐标迁移页把旧报告里的 SNV 坐标迁到新版参考序列。"""
from __future__ import annotations

import json
import os
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import Depends, FastAPI, File, Form, HTTPException, Request, UploadFile
from fastapi.responses import HTMLResponse, RedirectResponse
from fastapi.templating import Jinja2Templates

from . import db
from .compare import CompareError, categorize, parse_fasta, parse_vcf
from .migrate import MigrateError, MigratedRow, migrate_site, parse_mapping_tsv

TEMPLATES_DIR = Path(__file__).parent / "templates"
DEFAULT_DB_PATH = "./data/vcf_compare.db"

CATEGORY_TITLES = (("common", "共有位点"), ("a_only", "甲独有位点"), ("b_only", "乙独有位点"))


def create_app(db_path: str | None = None) -> FastAPI:
    path = db_path or os.environ.get("DATABASE_PATH", DEFAULT_DB_PATH)

    @asynccontextmanager
    async def lifespan(app: FastAPI):
        db.init_db(app.state.db_path)
        yield

    app = FastAPI(title="肿瘤变异位点比对", lifespan=lifespan)
    app.state.db_path = path
    db.init_db(path)  # 建表幂等，保证测试/直接运行时库已就绪
    templates = Jinja2Templates(directory=str(TEMPLATES_DIR))

    def get_conn():
        conn = db.connect(app.state.db_path)
        try:
            yield conn
        finally:
            conn.close()

    @app.get("/", response_class=HTMLResponse)
    def index(request: Request, conn=Depends(get_conn)):
        return templates.TemplateResponse(
            request, "index.html", {"reports": db.list_reports(conn)}
        )

    @app.post("/reports")
    def create_report(
        fasta: UploadFile = File(...),
        vcf_a: UploadFile = File(...),
        vcf_b: UploadFile = File(...),
        conn=Depends(get_conn),
    ):
        meta = {
            "fasta_filename": fasta.filename,
            "vcf_a_filename": vcf_a.filename,
            "vcf_b_filename": vcf_b.filename,
            "ref_name": None,
            "ref_length": None,
            "ref_sha256": None,
        }
        try:
            # 三个文件的全部内容只进入本次请求的局部变量
            ref_text = _read_text(fasta, "参考序列 FASTA")
            a_text = _read_text(vcf_a, "甲 VCF")
            b_text = _read_text(vcf_b, "乙 VCF")
            ref_name, ref_seq = parse_fasta(ref_text)
            meta.update(
                ref_name=ref_name,
                ref_length=len(ref_seq),
                ref_sha256=db.ref_digest(ref_seq),
            )
            sites_a = parse_vcf(a_text, ref_seq, ref_name, "甲")
            sites_b = parse_vcf(b_text, ref_seq, ref_name, "乙")
            compared = categorize(sites_a, sites_b)
        except CompareError as exc:
            report_id = db.insert_report(
                conn, status="failed", error=str(exc), compared=[], **meta
            )
        else:
            report_id = db.insert_report(
                conn, status="done", error=None, compared=compared, **meta
            )
        return RedirectResponse(url=f"/reports/{report_id}", status_code=303)

    @app.get("/reports/{report_id}", response_class=HTMLResponse)
    def report_detail(report_id: int, request: Request, conn=Depends(get_conn)):
        report, variants = db.get_report(conn, report_id)
        if report is None:
            raise HTTPException(status_code=404, detail="报告不存在")
        sections = [
            {
                "key": key,
                "title": title,
                "rows": [_variant_view(v) for v in variants if v["category"] == key],
            }
            for key, title in CATEGORY_TITLES
        ]
        return templates.TemplateResponse(
            request, "report.html", {"report": report, "sections": sections}
        )

    def _migrate_context(error: str | None, conn) -> dict:
        return {
            "reports": db.list_migratable_reports(conn),
            "migrations": db.list_migrations(conn),
            "error": error,
        }

    @app.get("/migrate", response_class=HTMLResponse)
    def migrate_page(request: Request, conn=Depends(get_conn)):
        return templates.TemplateResponse(
            request, "migrate.html", _migrate_context(None, conn)
        )

    @app.post("/migrate", response_class=HTMLResponse)
    def create_migration(
        request: Request,
        source_report_id: int = Form(...),
        old_fasta: UploadFile = File(...),
        new_fasta: UploadFile = File(...),
        mapping: UploadFile = File(...),
        conn=Depends(get_conn),
    ):
        report, variants = db.get_report(conn, source_report_id)
        if report is None:
            raise HTTPException(status_code=404, detail="源报告不存在")
        try:
            # 源报告须已完成、含 SNV 且无插入缺失
            snv_rows = [
                v for v in variants if len(v["ref"]) == 1 and len(v["alt"]) == 1
            ]
            if (
                report["status"] != "done"
                or not snv_rows
                or len(snv_rows) != len(variants)
            ):
                raise MigrateError(
                    f"报告 #{source_report_id} 不能迁移："
                    "源报告须已完成、含 SNV 且无插入缺失"
                )
            # 三个文件的全部内容只进入本次请求的局部变量
            old_text = _read_text(old_fasta, "旧参考序列 FASTA")
            new_text = _read_text(new_fasta, "新参考序列 FASTA")
            map_text = _read_text(mapping, "映射 TSV")
            old_name, old_seq = parse_fasta(old_text)
            new_name, new_seq = parse_fasta(new_text)
            old_sha256 = db.ref_digest(old_seq)
            if old_sha256 != report["ref_sha256"]:
                raise MigrateError(
                    f"旧参考序列 SHA-256（{old_sha256[:16]}…）与报告 "
                    f"#{source_report_id} 记录的参考摘要"
                    f"（{report['ref_sha256'][:16]}…）不一致"
                )
            block_map = parse_mapping_tsv(map_text, len(old_seq), len(new_seq))
            rows = []
            for v in snv_rows:
                site = migrate_site(
                    block_map, new_seq, new_name, v["pos"], v["ref"], v["alt"]
                )
                rows.append(
                    MigratedRow(
                        site=site,
                        category=v["category"],
                        lines_a=tuple(json.loads(v["lines_a"])),
                        lines_b=tuple(json.loads(v["lines_b"])),
                    )
                )
        except (CompareError, MigrateError) as exc:
            # 校验/迁移失败：不生成任何迁移报告，仅在本页提示原因
            return templates.TemplateResponse(
                request,
                "migrate.html",
                _migrate_context(str(exc), conn),
                status_code=400,
            )
        migration_id, _created = db.insert_migration(
            conn,
            source_report_id=source_report_id,
            old_ref_name=old_name,
            old_ref_length=len(old_seq),
            old_ref_sha256=old_sha256,
            new_ref_name=new_name,
            new_ref_length=len(new_seq),
            new_ref_sha256=db.ref_digest(new_seq),
            n_blocks=len(block_map),
            mapping_sha256=block_map.digest(),
            rows=rows,
        )
        return RedirectResponse(url=f"/migrations/{migration_id}", status_code=303)

    @app.get("/migrations/{migration_id}", response_class=HTMLResponse)
    def migration_detail(migration_id: int, request: Request, conn=Depends(get_conn)):
        migration, variants = db.get_migration(conn, migration_id)
        if migration is None:
            raise HTTPException(status_code=404, detail="迁移报告不存在")
        sections = [
            {
                "key": key,
                "title": title,
                "rows": [_migrated_view(v) for v in variants if v["category"] == key],
            }
            for key, title in CATEGORY_TITLES
        ]
        return templates.TemplateResponse(
            request,
            "migration.html",
            {"migration": migration, "sections": sections},
        )

    return app


def _read_text(upload: UploadFile, label: str) -> str:
    try:
        return upload.file.read().decode("utf-8-sig")
    except UnicodeDecodeError:
        raise CompareError(f"{label} 文件不是有效的 UTF-8 文本")


def _variant_view(row) -> dict:
    return {
        "pos": row["pos"],
        "ref": row["ref"] or "-",
        "alt": row["alt"] or "-",
        "kind": _kind_label(row["ref"], row["alt"]),
        "lines_a": ", ".join(str(n) for n in json.loads(row["lines_a"])) or "-",
        "lines_b": ", ".join(str(n) for n in json.loads(row["lines_b"])) or "-",
    }


def _migrated_view(row) -> dict:
    return {
        "pos": row["pos"],
        "ref": row["ref"],
        "alt": row["alt"],
        "src_pos": row["src_pos"],
        "lines_a": ", ".join(str(n) for n in json.loads(row["lines_a"])) or "-",
        "lines_b": ", ".join(str(n) for n in json.loads(row["lines_b"])) or "-",
    }


def _kind_label(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    return "插入" if ref == "" else "缺失"


app = create_app()
