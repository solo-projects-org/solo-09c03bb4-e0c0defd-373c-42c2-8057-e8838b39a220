"""核心比对逻辑单元测试。"""
import pytest

from app.compare import (
    MAX_RECORDS,
    CompareError,
    Site,
    _normalize,
    categorize,
    parse_fasta,
    parse_vcf,
)

# 1-based:
# 1-10  GTCAGTACGT
# 11-20 AAAAAACGTA   （11-16 为 poly-A）
# 21-30 CGATCGTACG
# 31-40 TAGCTAGCTA
# 41-50 GCTAGCTAGC
# 51-60 ATCGATCGAT
# 61-64 GGCC
REF = "GTCAGTACGTAAAAAACGTACGATCGTACGTAGCTAGCTAGCTAGCTAGCATCGATCGATGGCC"
NAME = "chrT"
HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def vcf_line(pos, ref, alt, chrom=NAME):
    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t."


def parse(lines, ref=REF, name=NAME):
    return parse_vcf(HEADER + "\n".join(lines) + "\n", ref, name, "甲")


# ---------- FASTA ----------

def test_fasta_ok_multiline_lowercase():
    name, seq = parse_fasta(">chr1 some description\nACGT\nacgt nn\n")
    assert name == "chr1"
    assert seq == "ACGTACGTNN"


def test_fasta_rejects_multiple_sequences():
    with pytest.raises(CompareError):
        parse_fasta(">chr1\nACGT\n>chr2\nACGT\n")


def test_fasta_rejects_bad_char():
    with pytest.raises(CompareError):
        parse_fasta(">chr1\nACGX\n")


def test_fasta_rejects_missing_header():
    with pytest.raises(CompareError):
        parse_fasta("ACGT\n")


def test_fasta_rejects_empty_sequence():
    with pytest.raises(CompareError):
        parse_fasta(">chr1\n")


# ---------- 归一化：SNV ----------

def test_snv_passthrough():
    assert parse([vcf_line(4, "A", "G")]) == {Site(4, "A", "G"): [3]}


def test_snv_lowercase_normalized():
    assert parse([vcf_line(4, "a", "g")]) == {Site(4, "A", "G"): [3]}


def test_snv_with_padded_prefix():
    # REF=CA ALT=CG（pos 3：ref[2:4]="CA"）等价于 pos 4 A>G
    assert parse([vcf_line(3, "CA", "CG")]) == {Site(4, "A", "G"): [3]}


# ---------- 归一化：插入缺失 ----------

def test_insertion_left_aligned_in_homopolymer():
    # pos 16 处插 A，位于 11-16 的 poly-A 中，应左对齐到 pos 11 之前
    assert parse([vcf_line(16, "A", "AA")]) == {Site(11, "", "A"): [3]}


def test_deletion_left_aligned_in_homopolymer():
    # pos 15-16 删一个 A，应左对齐到 pos 11
    assert parse([vcf_line(15, "AA", "A")]) == {Site(11, "A", ""): [3]}


def test_equivalent_writings_same_site():
    # 同一插入的两种写法（pos 13 与 pos 16 各插一个 A）归并为同一位点
    sites = parse([vcf_line(13, "A", "AA"), vcf_line(16, "A", "AA")])
    assert sites == {Site(11, "", "A"): [3, 4]}


def test_suffix_trim_insertion():
    # REF=GT ALT=GTT（pos 30：ref[29:31]="GT"）等价于 pos 31 前插 T
    assert parse([vcf_line(30, "GT", "GTT")]) == {Site(31, "", "T"): [3]}


def test_prefix_trim_deletion():
    # REF=AT ALT=A（pos 51：ref[50:52]="AT"）等价于 pos 52 删 T
    assert parse([vcf_line(51, "AT", "A")]) == {Site(52, "T", ""): [3]}


def test_normalize_deletion_rotation():
    # 参考序列 GACACAC（AC 串联重复）：删 pos 4-5 的 AC 与删 pos 2-3 的 AC
    # 产生相同突变序列，应左对齐到同一代表
    ref_seq = "GACACAC"
    assert _normalize(ref_seq, 4, "ACA", "A") == Site(2, "AC", "")
    assert _normalize(ref_seq, 2, "ACA", "A") == Site(2, "AC", "")
    # GACAC 中删 pos 5 的 C 与删 pos 3 的 C 突变结果不同，不得左移混同
    assert _normalize("GACAC", 4, "AC", "A") == Site(5, "C", "")


def test_normalize_insertion_rotation():
    # 同一参考序列上的插入同样左对齐到串联重复最左端
    ref_seq = "GACACAC"
    assert _normalize(ref_seq, 5, "C", "CAC") == Site(2, "", "AC")


# ---------- 非法记录 ----------

def test_ref_mismatch_rejected():
    with pytest.raises(CompareError, match="与参考序列不符"):
        parse([vcf_line(4, "C", "G")])


def test_pos_out_of_range_rejected():
    with pytest.raises(CompareError, match="超出参考序列范围"):
        parse([vcf_line(65, "C", "A")])


def test_pos_zero_rejected():
    with pytest.raises(CompareError, match="正整数"):
        parse([vcf_line(0, "A", "G")])


def test_multiallelic_rejected():
    with pytest.raises(CompareError, match="多等位"):
        parse([vcf_line(4, "A", "G,T")])


def test_symbolic_alt_rejected():
    with pytest.raises(CompareError):
        parse([vcf_line(4, "A", "<DEL>")])


def test_star_alt_rejected():
    with pytest.raises(CompareError):
        parse([vcf_line(4, "A", "*")])


def test_n_in_allele_rejected():
    with pytest.raises(CompareError):
        parse([vcf_line(4, "N", "G")])


def test_ref_eq_alt_rejected():
    with pytest.raises(CompareError, match="相同"):
        parse([vcf_line(4, "A", "A")])


def test_mnp_rejected():
    # pos 21-22 参考为 CG，换成 AT 是 MNP，不是 SNV/纯插入缺失
    with pytest.raises(CompareError, match="SNV 或纯插入缺失"):
        parse([vcf_line(21, "CG", "AT")])


def test_chrom_mismatch_rejected():
    with pytest.raises(CompareError, match="不一致"):
        parse([vcf_line(4, "A", "G", chrom="chrX")])


def test_too_many_records_rejected():
    lines = [vcf_line(4, "A", "G")] * (MAX_RECORDS + 1)
    with pytest.raises(CompareError, match="超过上限"):
        parse(lines)


def test_short_line_rejected():
    with pytest.raises(CompareError, match="列数不足"):
        parse(["chrT\t4\t.\tA"])


def test_error_message_contains_line_number():
    with pytest.raises(CompareError, match="第 4 行"):
        parse([vcf_line(4, "A", "G"), vcf_line(4, "C", "G")])


# ---------- 归类 ----------

def test_categorize():
    a = {Site(1, "A", "T"): [3], Site(5, "", "C"): [4, 7]}
    b = {Site(1, "A", "T"): [2], Site(9, "G", ""): [6]}
    rows = categorize(a, b)
    got = {(r.site, r.category, r.lines_a, r.lines_b) for r in rows}
    assert got == {
        (Site(1, "A", "T"), "common", (3,), (2,)),
        (Site(5, "", "C"), "a_only", (4, 7), ()),
        (Site(9, "G", ""), "b_only", (), (6,)),
    }
