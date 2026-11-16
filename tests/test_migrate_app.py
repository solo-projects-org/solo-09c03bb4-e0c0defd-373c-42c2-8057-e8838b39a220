"""坐标迁移网页集成测试：迁移流程、失败不生成报告、幂等去重、并发。"""
import asyncio
import re
from pathlib import Path

import httpx
import pytest
from fastapi.testclient import TestClient

from app.main import create_app

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"
REF_V1 = (EXAMPLES / "ref.fa").read_bytes()
REF_V2 = (EXAMPLES / "ref_v2.fa").read_bytes()
MAP_TSV = (EXAMPLES / "map_v1_v2.tsv").read_bytes()
A_SNV = (EXAMPLES / "lab_a_snv.vcf").read_bytes()
B_SNV = (EXAMPLES / "lab_b_snv.vcf").read_bytes()


@pytest.fixture()
def client(tmp_path):
    app = create_app(str(tmp_path / "test.db"))
    with TestClient(app) as c:
        yield c


def _create_report(client, a=A_SNV, b=B_SNV, fasta=REF_V1) -> int:
    resp = client.post(
        "/reports",
        files={
            "fasta": ("ref.fa", fasta, "text/plain"),
            "vcf_a": ("a.vcf", a, "text/plain"),
            "vcf_b": ("b.vcf", b, "text/plain"),
        },
        follow_redirects=False,
    )
    assert resp.status_code == 303
    return int(resp.headers["location"].rsplit("/", 1)[1])


def _migrate(client, report_id, old=REF_V1, new=REF_V2, mapping=MAP_TSV):
    return client.post(
        "/migrate",
        data={"source_report_id": str(report_id)},
        files={
            "old_fasta": ("old.fa", old, "text/plain"),
            "new_fasta": ("new.fa", new, "text/plain"),
            "mapping": ("map.tsv", mapping, "text/plain"),
        },
        follow_redirects=False,
    )


def _snv_report(client) -> int:
    """建一份仅含 SNV 的已完成报告（共有 2 / 甲独有 1 / 乙独有 1）。"""
    return _create_report(client)


# ---------- 页面与成功流程 ----------

def test_migrate_page_lists_only_snv_only_done_reports(client):
    snv_id = _snv_report(client)
    indel_id = _create_report(
        client,
        a=(EXAMPLES / "lab_a.vcf").read_bytes(),  # 含插入缺失
        b=(EXAMPLES / "lab_b.vcf").read_bytes(),
    )
    page = client.get("/migrate")
    assert page.status_code == 200
    assert f'value="{snv_id}"' in page.text
    assert f'value="{indel_id}"' not in page.text


def test_successful_migration(client):
    rid = _snv_report(client)
    resp = _migrate(client, rid)
    assert resp.status_code == 303
    assert resp.headers["location"] == "/migrations/1"

    page = client.get("/migrations/1")
    assert page.status_code == 200
    assert "迁移报告 #1" in page.text
    assert f'href="/reports/{rid}"' in page.text  # 标明源报告
    assert "新参考序列" in page.text and "sha256" in page.text  # 新参考摘要
    # 归类保留：共有 2 / 甲独有 1 / 乙独有 1
    assert "共有位点（2）" in page.text
    assert "甲独有位点（1）" in page.text
    assert "乙独有位点（1）" in page.text
    # 正向块：4 A>G → 4 A>G（甲行 4 / 乙行 3）
    assert re.search(
        r"<td>4</td>\s*<td>A</td>\s*<td>G</td>\s*<td>4</td>\s*<td>4</td>\s*<td>3</td>",
        page.text,
    )
    # 负向块：40 A>T → 21 T>A（反向互补，甲行 6 / 乙行 5）
    assert re.search(
        r"<td>21</td>\s*<td>T</td>\s*<td>A</td>\s*<td>40</td>\s*<td>6</td>\s*<td>5</td>",
        page.text,
    )
    # 甲独有 20 A>C → 20 A>C（甲行 5）
    assert re.search(r"<td>20</td>\s*<td>A</td>\s*<td>C</td>\s*<td>20</td>\s*<td>5</td>", page.text)
    # 乙独有 25 C>T → 36 G>A（反向互补，乙行 4）
    assert re.search(r"<td>36</td>\s*<td>G</td>\s*<td>A</td>\s*<td>25</td>\s*<td>4</td>", page.text)
    # 迁移页历史列表可见
    assert 'href="/migrations/1"' in client.get("/migrate").text


def test_old_report_untouched_after_migration(client):
    rid = _snv_report(client)
    _migrate(client, rid)
    page = client.get(f"/reports/{rid}")
    assert page.status_code == 200
    assert "共有位点（2）" in page.text
    assert "<td>40</td>" in page.text  # 旧坐标不变


# ---------- 失败：不生成迁移报告 ----------

def _assert_no_migration(client):
    assert "还没有迁移报告" in client.get("/migrate").text


def test_old_fasta_sha_mismatch(client):
    rid = _snv_report(client)
    resp = _migrate(client, rid, old=REF_V2)  # 用错旧参考
    assert resp.status_code == 400
    assert "SHA-256" in resp.text and "不一致" in resp.text
    _assert_no_migration(client)


def test_unmapped_site_fails(client):
    rid = _snv_report(client)
    mapping = b"src_start\tsrc_end\tdst_start\tdst_end\tstrand\n1\t20\t1\t20\t+\n"
    resp = _migrate(client, rid, mapping=mapping)
    assert resp.status_code == 400
    assert "不在任何映射块" in resp.text
    assert "25" in resp.text  # 第一个无映射位点
    _assert_no_migration(client)


def test_target_ref_mismatch_fails(client):
    rid = _snv_report(client)
    # 把新参考第 36 位（G）改成 A，使 25 C>T 迁移后的 REF 核对失败
    seq = bytearray(REF_V2.decode().split("\n")[1], "ascii")
    assert seq[35:36] == b"G"
    seq[35:36] = b"A"
    bad_new = f">chrT\n{seq.decode()}\n".encode()
    resp = _migrate(client, rid, new=bad_new)
    assert resp.status_code == 400
    assert "与新参考序列不符" in resp.text
    _assert_no_migration(client)


def test_bad_mapping_tsv_fails(client):
    rid = _snv_report(client)
    resp = _migrate(client, rid, mapping=b"1\t20\t1\t19\t+\n")  # 不等长
    assert resp.status_code == 400
    assert "长度不等" in resp.text
    _assert_no_migration(client)


def test_report_with_indels_rejected(client):
    rid = _create_report(
        client,
        a=(EXAMPLES / "lab_a.vcf").read_bytes(),
        b=(EXAMPLES / "lab_b.vcf").read_bytes(),
    )
    resp = _migrate(client, rid)
    assert resp.status_code == 400
    assert "不能迁移" in resp.text
    _assert_no_migration(client)


def test_failed_report_rejected(client):
    rid = _create_report(client, a=(EXAMPLES / "lab_bad.vcf").read_bytes())
    resp = _migrate(client, rid)
    assert resp.status_code == 400
    assert "不能迁移" in resp.text
    _assert_no_migration(client)


def test_missing_source_report(client):
    resp = _migrate(client, 999)
    assert resp.status_code == 404


def test_migration_not_found(client):
    assert client.get("/migrations/999").status_code == 404


# ---------- 幂等与并发 ----------

def test_duplicate_submission_returns_same_migration(client):
    rid = _snv_report(client)
    r1 = _migrate(client, rid)
    r2 = _migrate(client, rid)
    assert r1.headers["location"] == "/migrations/1"
    assert r2.headers["location"] == "/migrations/1"
    # 历史列表只有一条迁移记录
    assert client.get("/migrate").text.count('href="/migrations/') == 1


def test_different_mapping_content_creates_new_migration(client):
    rid = _snv_report(client)
    _migrate(client, rid)
    # 同样的块拆成更细的等长块，内容不同 → 新迁移报告
    mapping = (
        b"src_start\tsrc_end\tdst_start\tdst_end\tstrand\n"
        b"1\t10\t1\t10\t+\n11\t20\t11\t20\t+\n"
        b"21\t40\t21\t40\t-\n41\t64\t41\t64\t+\n"
    )
    resp = _migrate(client, rid, mapping=mapping)
    assert resp.headers["location"] == "/migrations/2"


def test_concurrent_identical_submissions_dedup(tmp_path):
    """4 个完全相同的迁移请求并发提交：只产生一份完整迁移报告。"""
    app = create_app(str(tmp_path / "concurrent.db"))
    rid = None

    async def run():
        transport = httpx.ASGITransport(app=app)
        async with httpx.AsyncClient(
            transport=transport, base_url="http://t", follow_redirects=False
        ) as ac:
            nonlocal rid
            resp = await ac.post(
                "/reports",
                files={
                    "fasta": ("ref.fa", REF_V1, "text/plain"),
                    "vcf_a": ("a.vcf", A_SNV, "text/plain"),
                    "vcf_b": ("b.vcf", B_SNV, "text/plain"),
                },
            )
            rid = resp.headers["location"].rsplit("/", 1)[1]

            async def post():
                return await ac.post(
                    "/migrate",
                    data={"source_report_id": rid},
                    files={
                        "old_fasta": ("old.fa", REF_V1, "text/plain"),
                        "new_fasta": ("new.fa", REF_V2, "text/plain"),
                        "mapping": ("map.tsv", MAP_TSV, "text/plain"),
                    },
                )

            return await asyncio.gather(*(post() for _ in range(4)))

    resps = asyncio.run(run())
    locations = {r.headers["location"] for r in resps}
    assert locations == {"/migrations/1"}  # 全部指向同一份迁移报告

    with TestClient(app) as c:
        page = c.get("/migrate")
        assert page.text.count('href="/migrations/') == 1  # 无重复
        detail = c.get("/migrations/1")
        # 无半份报告：4 个位点全部迁移，归类完整
        assert "共有位点（2）" in detail.text
        assert "甲独有位点（1）" in detail.text
        assert "乙独有位点（1）" in detail.text
