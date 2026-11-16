"""坐标迁移核心逻辑单元测试。"""
import pytest

from app.migrate import (
    MAX_BLOCKS,
    MigrateError,
    Block,
    BlockMap,
    MigratedSite,
    migrate_site,
    parse_mapping_tsv,
    reverse_complement,
)

# 1-based:
# 1-20  GTCAGTACGTAAAAAACGTA
# 21-40 TAGCTAGCTACGTACGATCG   （旧 21-40 的反向互补）
# 41-64 GCTAGCTAGCATCGATCGATGGCC
NEW = "GTCAGTACGTAAAAAACGTATAGCTAGCTACGTACGATCGGCTAGCTAGCATCGATCGATGGCC"
NEW_NAME = "chrT"
MAP_TEXT = (
    "src_start\tsrc_end\tdst_start\tdst_end\tstrand\n"
    "1\t20\t1\t20\t+\n"
    "21\t40\t21\t40\t-\n"
    "41\t64\t41\t64\t+\n"
)


def parse(text=MAP_TEXT, src_len=64, dst_len=64):
    return parse_mapping_tsv(text, src_len, dst_len)


# ---------- 反向互补 ----------

def test_reverse_complement():
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("A") == "T"
    assert reverse_complement("C") == "G"
    assert reverse_complement("AAAC") == "GTTT"


# ---------- 映射 TSV 解析 ----------

def test_parse_ok_with_header_and_comments():
    bm = parse("# 注释行\n\n" + MAP_TEXT)
    assert len(bm) == 3
    assert bm.blocks[0] == Block(1, 20, 1, 20, "+")
    assert bm.blocks[1] == Block(21, 40, 21, 40, "-")
    assert bm.blocks[2] == Block(41, 64, 41, 64, "+")


def test_parse_ok_without_header():
    bm = parse("1\t20\t1\t20\t+\n21\t40\t21\t40\t-\n41\t64\t41\t64\t+\n")
    assert len(bm) == 3


def test_digest_stable_across_row_order_and_header():
    a = parse(MAP_TEXT)
    shuffled = "41\t64\t41\t64\t+\n1\t20\t1\t20\t+\n21\t40\t21\t40\t-\n"
    b = parse(shuffled)
    assert a.digest() == b.digest()


def test_digest_changes_with_content():
    a = parse(MAP_TEXT)
    b = parse("1\t20\t1\t20\t+\n21\t40\t21\t40\t-\n41\t64\t41\t64\t+\n".replace(
        "41\t64\t41\t64", "41\t63\t41\t63"))
    # 上面替换后第三块不再覆盖 64，但仍是合法映射，摘要应不同
    assert a.digest() != b.digest()


def test_reject_empty():
    with pytest.raises(MigrateError, match="没有任何映射块"):
        parse("src_start\tsrc_end\tdst_start\tdst_end\tstrand\n")


def test_reject_bad_column_count():
    with pytest.raises(MigrateError, match="第 2 行.*5 列"):
        parse("1\t20\t1\t20\t+\n2\t3\t4\n")


def test_reject_non_integer():
    with pytest.raises(MigrateError, match="第 2 行.*正整数"):
        parse("1\t20\t1\t20\t+\n2\tx\t2\t3\t+\n")


def test_reject_inverted_interval():
    with pytest.raises(MigrateError, match="起止颠倒"):
        parse("20\t1\t1\t20\t+\n")


def test_reject_unequal_length():
    with pytest.raises(MigrateError, match="长度不等"):
        parse("1\t20\t1\t19\t+\n")


def test_reject_bad_strand():
    with pytest.raises(MigrateError, match="方向"):
        parse("1\t20\t1\t20\t?\n")


def test_reject_src_out_of_range():
    with pytest.raises(MigrateError, match="旧参考序列范围"):
        parse("1\t65\t1\t65\t+\n")


def test_reject_dst_out_of_range():
    with pytest.raises(MigrateError, match="新参考序列范围"):
        parse("1\t64\t2\t65\t+\n")


def test_reject_src_overlap():
    with pytest.raises(MigrateError, match="源区间重叠"):
        parse("1\t20\t1\t20\t+\n10\t30\t31\t51\t+\n")


def test_reject_dst_overlap():
    with pytest.raises(MigrateError, match="目标区间重叠"):
        parse("1\t20\t1\t20\t+\n21\t40\t20\t39\t+\n")


def test_reject_too_many_blocks():
    text = "".join(f"{i}\t{i}\t{i}\t{i}\t+\n" for i in range(1, MAX_BLOCKS + 2))
    with pytest.raises(MigrateError, match="超过上限"):
        parse(text, src_len=MAX_BLOCKS + 1, dst_len=MAX_BLOCKS + 1)


# ---------- 逐位点映射 ----------

def migrate(pos, ref, alt, bm=None, new_seq=NEW):
    return migrate_site(bm or parse(), new_seq, NEW_NAME, pos, ref, alt)


def test_plus_strand_passthrough():
    assert migrate(4, "A", "G") == MigratedSite(src_pos=4, pos=4, ref="A", alt="G")
    assert migrate(20, "A", "C") == MigratedSite(src_pos=20, pos=20, ref="A", alt="C")


def test_minus_strand_reverse_complements():
    # 旧 25 C>T → 新 36，C/T 互补为 G/A
    assert migrate(25, "C", "T") == MigratedSite(src_pos=25, pos=36, ref="G", alt="A")
    # 旧 40 A>T → 新 21，A/T 互补为 T/A
    assert migrate(40, "A", "T") == MigratedSite(src_pos=40, pos=21, ref="T", alt="A")


def test_block_edges():
    # 负向块首端映射到目标区间末端，末端映射到首端
    # 源 21（C）→ 目标 40，互补为 G，新参考 40 位正是 G
    assert migrate(21, "C", "A") == MigratedSite(src_pos=21, pos=40, ref="G", alt="T")
    # 源 40（A）→ 目标 21，互补为 T
    assert migrate(40, "A", "T").pos == 21


def test_unmapped_site_rejected():
    bm = parse("1\t20\t1\t20\t+\n")
    with pytest.raises(MigrateError, match="不在任何映射块"):
        migrate(30, "C", "T", bm=bm)


def test_target_ref_mismatch_rejected():
    # 新参考 36 位是 G，若声称源 25 位 REF 为 A（互补 T）则核对失败
    with pytest.raises(MigrateError, match="与新参考序列不符"):
        migrate(25, "A", "T")


def test_mapping_is_unique_per_site():
    # 源区间互不重叠 ⇒ 任一坐标至多一个覆盖块
    bm = parse()
    assert bm.find(20).src_start == 1
    assert bm.find(21).src_start == 21
    assert bm.find(64).src_start == 41
    assert BlockMap([]).find(1) is None
