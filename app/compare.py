"""核心比对逻辑：FASTA/VCF 解析、REF 核对、插入缺失左对齐归一化、位点集合比较。

本模块不做任何 I/O，全部基于传入的文本与参考序列字符串，保证一次上传的
参考序列、位点、行号只存在于本次调用的局部变量中，从机制上避免并发批次混用。
"""
from __future__ import annotations

from dataclasses import dataclass

BASES = frozenset("ACGT")
REF_BASES = frozenset("ACGTN")
MAX_RECORDS = 1000  # 每份 VCF 最多 1000 条记录


class CompareError(Exception):
    """任一非法记录 / 非法参考序列导致的整次比对失败。"""


@dataclass(frozen=True, order=True)
class Site:
    """归一化后的位点（等价写法的唯一代表）。

    pos: 1-based，去除共有前后缀并左对齐后的起始位置。
    ref: 修剪后的参考等位，纯插入时为 ""。
    alt: 修剪后的变异等位，纯缺失时为 ""。
    """

    pos: int
    ref: str
    alt: str

    @property
    def kind(self) -> str:
        if len(self.ref) == 1 and len(self.alt) == 1:
            return "SNV"
        if self.ref == "":
            return "插入"
        return "缺失"


@dataclass(frozen=True)
class ComparedSite:
    """一个位点在两份 VCF 中的归属与来源行号。"""

    site: Site
    category: str  # "common" | "a_only" | "b_only"
    lines_a: tuple[int, ...]
    lines_b: tuple[int, ...]


def parse_fasta(text: str) -> tuple[str, str]:
    """解析单条染色体的 FASTA，返回 (序列名, 大写序列)。非法则抛 CompareError。"""
    name: str | None = None
    chunks: list[str] = []
    for lineno, raw in enumerate(text.split("\n"), 1):
        line = raw.rstrip("\r").strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                raise CompareError(
                    f"参考序列 FASTA 第 {lineno} 行出现第二条序列，一次只能上传单条染色体"
                )
            header = line[1:].strip()
            if not header:
                raise CompareError(f"参考序列 FASTA 第 {lineno} 行缺少序列名")
            name = header.split()[0]
            continue
        if name is None:
            raise CompareError("参考序列 FASTA 首行必须以 '>' 开头")
        seq = line.replace(" ", "").replace("\t", "").upper()
        bad = sorted(set(seq) - REF_BASES)
        if bad:
            raise CompareError(
                f"参考序列 FASTA 第 {lineno} 行含非法字符: {''.join(bad)}"
            )
        chunks.append(seq)
    if name is None:
        raise CompareError("参考序列 FASTA 为空或缺少 '>' 头行")
    sequence = "".join(chunks)
    if not sequence:
        raise CompareError("参考序列 FASTA 没有任何序列内容")
    return name, sequence


def _normalize(ref_seq: str, pos: int, ref: str, alt: str) -> Site | None:
    """去除共有前后缀并按参考序列左对齐，返回归一化位点。

    仅接受 SNV 与纯插入缺失；其余（如 MNP、复杂替换）返回 None。
    pos 为 1-based，ref/alt 为大写 ACGT 字符串，且已确认 ref 与参考序列一致。
    """
    # 去共有后缀
    while ref and alt and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]
    # 去共有前缀
    while ref and alt and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos += 1
    if len(ref) == 1 and len(alt) == 1:
        return Site(pos, ref, alt)  # SNV 无需左对齐
    if ref and alt:
        return None  # MNP / 复杂替换，不在允许范围
    # 纯插入（ref == ""）或纯缺失（alt == ""）：按参考序列向左滑动
    allele = ref if alt == "" else alt
    while pos > 1 and ref_seq[pos - 2] == allele[-1]:
        allele = allele[-1] + allele[:-1]
        pos -= 1
    if alt == "":
        return Site(pos, allele, "")
    return Site(pos, "", allele)


def parse_vcf(
    text: str, ref_seq: str, ref_name: str, label: str
) -> dict[Site, list[int]]:
    """解析一份 VCF，返回 {归一化位点: [原文件行号, ...]}（升序）。

    同一文件内的等价位点合并并保留全部原行号；任一记录非法即抛 CompareError，
    由调用方终止整次比对。
    """
    sites: dict[Site, list[int]] = {}
    n_records = 0
    for lineno, raw in enumerate(text.split("\n"), 1):
        line = raw.rstrip("\r")
        if not line.strip() or line.startswith("#"):
            continue
        n_records += 1
        if n_records > MAX_RECORDS:
            raise CompareError(
                f"{label} VCF 记录数超过上限 {MAX_RECORDS} 条（第 {lineno} 行）"
            )
        cols = line.split("\t") if "\t" in line else line.split()
        if len(cols) < 5:
            raise CompareError(
                f"{label} VCF 第 {lineno} 行列数不足，至少需要 CHROM/POS/ID/REF/ALT 5 列"
            )
        chrom, pos_s, _id, ref, alt = cols[0], cols[1], cols[2], cols[3].upper(), cols[4].upper()
        if chrom != ref_name:
            raise CompareError(
                f"{label} VCF 第 {lineno} 行染色体 '{chrom}' 与参考序列 '{ref_name}' 不一致"
            )
        if not pos_s.isdigit() or int(pos_s) < 1:
            raise CompareError(f"{label} VCF 第 {lineno} 行 POS 不是正整数: '{pos_s}'")
        pos = int(pos_s)
        if not ref or any(c not in BASES for c in ref):
            raise CompareError(
                f"{label} VCF 第 {lineno} 行 REF 含非法字符（仅允许 A/C/G/T）: '{cols[3]}'"
            )
        if "," in alt:
            raise CompareError(
                f"{label} VCF 第 {lineno} 行 ALT 为多等位 '{cols[4]}'，仅支持双等位"
            )
        if not alt or any(c not in BASES for c in alt):
            raise CompareError(
                f"{label} VCF 第 {lineno} 行 ALT 含非法字符（仅允许 A/C/G/T 的 SNV 或纯插入缺失）: '{cols[4]}'"
            )
        if ref == alt:
            raise CompareError(f"{label} VCF 第 {lineno} 行 REF 与 ALT 相同，不是变异")
        if pos + len(ref) - 1 > len(ref_seq):
            raise CompareError(
                f"{label} VCF 第 {lineno} 行 POS {pos} 加 REF 长度超出参考序列范围（长度 {len(ref_seq)}）"
            )
        actual = ref_seq[pos - 1 : pos - 1 + len(ref)]
        if actual != ref:
            raise CompareError(
                f"{label} VCF 第 {lineno} 行 REF '{ref}' 与参考序列不符（{ref_name}:{pos} 实际为 '{actual}'）"
            )
        site = _normalize(ref_seq, pos, ref, alt)
        if site is None:
            raise CompareError(
                f"{label} VCF 第 {lineno} 行不是 SNV 或纯插入缺失（REF={ref} ALT={alt}）"
            )
        sites.setdefault(site, []).append(lineno)
    return sites


def categorize(
    sites_a: dict[Site, list[int]], sites_b: dict[Site, list[int]]
) -> list[ComparedSite]:
    """按 共有 / 甲独有 / 乙独有 归类，按坐标排序。"""
    rows: list[ComparedSite] = []
    for site in sorted(set(sites_a) | set(sites_b)):
        la = tuple(sites_a.get(site, ()))
        lb = tuple(sites_b.get(site, ()))
        if la and lb:
            category = "common"
        elif la:
            category = "a_only"
        else:
            category = "b_only"
        rows.append(ComparedSite(site, category, la, lb))
    return rows
