"""坐标迁移：映射 TSV 解析校验、链向感知的逐位点映射、目标 REF 核对。

本模块不做任何 I/O：映射块、参考序列、位点全部来自调用方传入的字符串，
一次迁移的数据只存在于本次调用的局部变量中。

映射 TSV 每行 5 列（制表符或空白分隔）：

    src_start  src_end  dst_start  dst_end  strand

区间为 1 基闭区间；strand 为 '+' 或 '-'。允许空行、'#' 注释行；
若第一个有效行的前 4 列不全是正整数，则视为表头跳过。
"""
from __future__ import annotations

import bisect
import hashlib
from dataclasses import dataclass

MAX_BLOCKS = 10000  # 映射块数上限

_RC_TABLE = str.maketrans("ACGT", "TGCA")


class MigrateError(Exception):
    """映射文件非法或任一位点无法迁移（无映射 / 目标 REF 不符）导致的整体失败。"""


def reverse_complement(seq: str) -> str:
    """A/C/G/T 序列的反向互补。"""
    return seq.translate(_RC_TABLE)[::-1]


@dataclass(frozen=True)
class Block:
    """一个映射块：源区间与目标区间（均 1 基闭区间）等长，strand 为 '+' 或 '-'。"""

    src_start: int
    src_end: int
    dst_start: int
    dst_end: int
    strand: str


@dataclass(frozen=True)
class MigratedSite:
    """一个位点迁移后的结果：新坐标、（可能反向互补过的）REF/ALT 与原坐标。"""

    src_pos: int
    pos: int
    ref: str
    alt: str


@dataclass(frozen=True)
class MigratedRow:
    """迁移结果加上源报告的归类与甲乙 VCF 原行号，用于写入迁移报告。"""

    site: MigratedSite
    category: str  # "common" | "a_only" | "b_only"
    lines_a: tuple[int, ...]
    lines_b: tuple[int, ...]


class BlockMap:
    """按源区间排序、互不重叠的映射块集合，支持按源坐标查找唯一覆盖块。"""

    def __init__(self, blocks: list[Block]):
        self.blocks = tuple(blocks)  # 已按 src_start 排序
        self._starts = [b.src_start for b in self.blocks]

    def __len__(self) -> int:
        return len(self.blocks)

    def find(self, pos: int) -> Block | None:
        """返回覆盖源坐标 pos 的唯一块（源区间互不重叠，故至多一个）。"""
        i = bisect.bisect_right(self._starts, pos) - 1
        if i < 0:
            return None
        block = self.blocks[i]
        return block if pos <= block.src_end else None

    def digest(self) -> str:
        """映射内容摘要：规范化（排序后）块文本的 SHA-256。

        与行序、表头、注释无关——内容相同的映射文件得到相同摘要。
        """
        canonical = "\n".join(
            f"{b.src_start}-{b.src_end}:{b.dst_start}-{b.dst_end}:{b.strand}"
            for b in self.blocks
        )
        return hashlib.sha256(canonical.encode("ascii")).hexdigest()


def parse_mapping_tsv(text: str, src_len: int, dst_len: int) -> BlockMap:
    """解析并校验映射 TSV，返回 BlockMap。任一非法即抛 MigrateError。

    校验：恰好 5 列；正整数且 start<=end；源/目标区间等长、各自不越界；
    源区间两两不重叠、目标区间两两不重叠；strand 为 '+' 或 '-'。
    """
    blocks: list[Block] = []
    seen_data = False
    for lineno, raw in enumerate(text.split("\n"), 1):
        line = raw.rstrip("\r").strip()
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t") if "\t" in line else line.split()
        # 第一个有效行若前 4 列不全是正整数，视为表头
        if not seen_data and not blocks and any(
            not c.isdigit() or int(c) < 1 for c in cols[:4]
        ):
            seen_data = True
            continue
        seen_data = True
        if len(cols) != 5:
            raise MigrateError(
                f"映射 TSV 第 {lineno} 行需要恰好 5 列"
                "（src_start src_end dst_start dst_end strand），"
                f"实际为 {len(cols)} 列"
            )
        s1, s2, d1, d2, strand = cols
        for name, val in (("src_start", s1), ("src_end", s2), ("dst_start", d1), ("dst_end", d2)):
            if not val.isdigit() or int(val) < 1:
                raise MigrateError(
                    f"映射 TSV 第 {lineno} 行 {name} 不是正整数: '{val}'"
                )
        b = Block(int(s1), int(s2), int(d1), int(d2), strand)
        if b.src_end < b.src_start:
            raise MigrateError(
                f"映射 TSV 第 {lineno} 行源区间起止颠倒（{b.src_start} > {b.src_end}）"
            )
        if b.dst_end < b.dst_start:
            raise MigrateError(
                f"映射 TSV 第 {lineno} 行目标区间起止颠倒（{b.dst_start} > {b.dst_end}）"
            )
        if b.src_end - b.src_start != b.dst_end - b.dst_start:
            raise MigrateError(
                f"映射 TSV 第 {lineno} 行源区间与目标区间长度不等"
                f"（{b.src_end - b.src_start + 1} ≠ {b.dst_end - b.dst_start + 1}）"
            )
        if b.src_end > src_len:
            raise MigrateError(
                f"映射 TSV 第 {lineno} 行源区间超出旧参考序列范围（长度 {src_len}）"
            )
        if b.dst_end > dst_len:
            raise MigrateError(
                f"映射 TSV 第 {lineno} 行目标区间超出新参考序列范围（长度 {dst_len}）"
            )
        if strand not in ("+", "-"):
            raise MigrateError(
                f"映射 TSV 第 {lineno} 行方向列须为 '+' 或 '-'，实际为 '{strand}'"
            )
        blocks.append(b)
        if len(blocks) > MAX_BLOCKS:
            raise MigrateError(f"映射 TSV 块数超过上限 {MAX_BLOCKS} 个（第 {lineno} 行）")
    if not blocks:
        raise MigrateError("映射 TSV 没有任何映射块")
    return BlockMap(_check_non_overlapping(blocks))


def _check_non_overlapping(blocks: list[Block]) -> list[Block]:
    """校验源、目标两侧区间各自两两不重叠，返回按源区间排序的块列表。"""
    by_src = sorted(blocks, key=lambda b: (b.src_start, b.src_end))
    for prev, cur in zip(by_src, by_src[1:]):
        if cur.src_start <= prev.src_end:
            raise MigrateError(
                f"映射块源区间重叠：{prev.src_start}-{prev.src_end} 与 "
                f"{cur.src_start}-{cur.src_end}"
            )
    by_dst = sorted(blocks, key=lambda b: (b.dst_start, b.dst_end))
    for prev, cur in zip(by_dst, by_dst[1:]):
        if cur.dst_start <= prev.dst_end:
            raise MigrateError(
                f"映射块目标区间重叠：{prev.dst_start}-{prev.dst_end} 与 "
                f"{cur.dst_start}-{cur.dst_end}"
            )
    return by_src


def migrate_site(
    block_map: BlockMap, new_seq: str, new_name: str, pos: int, ref: str, alt: str
) -> MigratedSite:
    """把旧参考坐标 pos 的位点（SNV）迁移到新参考坐标并核对目标 REF。

    正向块：新坐标 = dst_start + (pos - src_start)，REF/ALT 不变；
    负向块：新坐标 = dst_end - (pos - src_start)，REF/ALT 取反向互补。
    无覆盖块或新参考该处碱基与迁移后 REF 不符，抛 MigrateError。
    """
    block = block_map.find(pos)
    if block is None:
        raise MigrateError(
            f"源位点 {pos}（{ref}>{alt}）不在任何映射块的源区间内，无法迁移"
        )
    offset = pos - block.src_start
    if block.strand == "+":
        new_pos = block.dst_start + offset
        new_ref, new_alt = ref, alt
    else:
        new_pos = block.dst_end - offset
        new_ref, new_alt = reverse_complement(ref), reverse_complement(alt)
    actual = new_seq[new_pos - 1 : new_pos - 1 + len(new_ref)]
    if actual != new_ref:
        raise MigrateError(
            f"源位点 {pos} 迁移到 {new_name}:{new_pos} 后 REF '{new_ref}' "
            f"与新参考序列不符（实际为 '{actual}'）"
        )
    return MigratedSite(src_pos=pos, pos=new_pos, ref=new_ref, alt=new_alt)
