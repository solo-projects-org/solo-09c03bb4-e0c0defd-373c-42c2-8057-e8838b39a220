"""网页端集成测试：上传、报告查看、失败原子性、并发批次隔离。"""
import asyncio
from pathlib import Path

import httpx
import pytest
from fastapi.testclient import TestClient

from app.main import create_app

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"


@pytest.fixture()
def client(tmp_path):
    app = create_app(str(tmp_path / "test.db"))
    with TestClient(app) as c:
        yield c


def _files(fasta: bytes, a: bytes, b: bytes):
    return {
        "fasta": ("ref.fa", fasta, "text/plain"),
        "vcf_a": ("a.vcf", a, "text/plain"),
        "vcf_b": ("b.vcf", b, "text/plain"),
    }


def _example_files():
    return _files(
        (EXAMPLES / "ref.fa").read_bytes(),
        (EXAMPLES / "lab_a.vcf").read_bytes(),
        (EXAMPLES / "lab_b.vcf").read_bytes(),
    )


def test_index_loads(client):
    resp = client.get("/")
    assert resp.status_code == 200
    assert "上传新比对" in resp.text


def test_successful_comparison(client):
    resp = client.post("/reports", files=_example_files(), follow_redirects=False)
    assert resp.status_code == 303
    location = resp.headers["location"]

    page = client.get(location)
    assert page.status_code == 200
    # 4 共有 / 2 甲独有 / 2 乙独有
    assert "共有位点（4）" in page.text
    assert "甲独有位点（2）" in page.text
    assert "乙独有位点（2）" in page.text
    # 甲文件内两种等价写法合并为一个位点，保留全部原行号 8 和 9
    assert "8, 9" in page.text
    # 报告列表里也能看到
    index = client.get("/")
    assert "完成" in index.text


def test_invalid_record_fails_and_old_report_survives(client):
    ok = client.post("/reports", files=_example_files(), follow_redirects=False)
    good_url = ok.headers["location"]

    bad = _files(
        (EXAMPLES / "ref.fa").read_bytes(),
        (EXAMPLES / "lab_bad.vcf").read_bytes(),  # 第 4 行 REF 与参考序列不符
        (EXAMPLES / "lab_b.vcf").read_bytes(),
    )
    resp = client.post("/reports", files=bad, follow_redirects=False)
    assert resp.status_code == 303
    page = client.get(resp.headers["location"])
    assert "失败" in page.text
    assert "与参考序列不符" in page.text
    assert "第 4 行" in page.text

    # 旧报告不受影响，仍可查看
    old = client.get(good_url)
    assert old.status_code == 200
    assert "共有位点（4）" in old.text


def test_bad_fasta_fails(client):
    files = _files(b"not-a-fasta\n", b"#CHROM\n", b"#CHROM\n")
    resp = client.post("/reports", files=files, follow_redirects=False)
    page = client.get(resp.headers["location"])
    assert "失败" in page.text
    assert "FASTA" in page.text


def test_report_not_found(client):
    assert client.get("/reports/999").status_code == 404


def _batch(i: int):
    """构造第 i 批数据：参考序列名、三个位点坐标都按批次唯一。"""
    name = f"chrB{i}"
    seq = "ACGT" * 25  # 100 bp
    alt_of = {"A": "C", "C": "G", "G": "T", "T": "A"}
    base = i * 10
    p_common, p_a, p_b = base + 2, base + 5, base + 8

    def snv(p):
        ref = seq[p - 1]
        return f"{name}\t{p}\t.\t{ref}\t{alt_of[ref]}\t.\t.\t."

    header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    fasta = f">{name}\n{seq}\n".encode()
    a = (header + snv(p_common) + "\n" + snv(p_a) + "\n").encode()
    b = (header + snv(p_common) + "\n" + snv(p_b) + "\n").encode()
    return fasta, a, b, (p_common, p_a, p_b)


def test_concurrent_uploads_do_not_mix_batches(tmp_path):
    """6 批不同参考序列/位点并发出传，逐份报告核对没有串批次。"""
    app = create_app(str(tmp_path / "concurrent.db"))
    n = 6

    async def run():
        transport = httpx.ASGITransport(app=app)
        async with httpx.AsyncClient(
            transport=transport, base_url="http://t", follow_redirects=True
        ) as ac:
            async def post(i):
                fasta, a, b, _ = _batch(i)
                resp = await ac.post("/reports", files=_files(fasta, a, b))
                assert resp.status_code == 200  # 跟随重定向到报告页
                return resp

            return await asyncio.gather(*(post(i) for i in range(n)))

    pages = asyncio.run(run())
    seen_urls = set()
    for i, page in enumerate(pages):
        _, _, _, (p_common, p_a, p_b) = _batch(i)
        assert f"chrB{i}" in page.text, f"报告参考序列名串批次: {page.url}"
        assert "共有位点（1）" in page.text
        assert "甲独有位点（1）" in page.text
        assert "乙独有位点（1）" in page.text
        for p in (p_common, p_a, p_b):
            assert f"<td>{p}</td>" in page.text
        # 其他批次的坐标不得出现
        for j in range(n):
            if j == i:
                continue
            for p in _batch(j)[3]:
                assert f"<td>{p}</td>" not in page.text
        seen_urls.add(str(page.url))
    assert len(seen_urls) == n  # 6 份独立报告
