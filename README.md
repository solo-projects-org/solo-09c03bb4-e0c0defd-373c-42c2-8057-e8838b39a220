# 肿瘤变异位点比对

肿瘤变异数据员用于核对两家实验室（甲、乙）检测结果差异的网页工具。
一次上传**单条染色体的 FASTA 参考序列**和**甲、乙两份 VCF**，系统核对 REF、
归一化位点后列出**共有 / 甲独有 / 乙独有**位点及各自来源行号。

## 运行

### Docker Compose（推荐）

```bash
docker compose up --build
```

打开 http://localhost:8000 。数据保存在名为 `vcf_data` 的卷里（`/data/vcf_compare.db`）。

### 本地开发

```bash
pip install -r requirements.txt
uvicorn app.main:app --reload
# 打开 http://127.0.0.1:8000
```

数据库路径可用环境变量 `DATABASE_PATH` 指定（默认 `./data/vcf_compare.db`）。

### 测试

```bash
pip install -r requirements-dev.txt
pytest tests/
```

`examples/` 目录有一组可直接上传的示例文件（`ref.fa` + `lab_a.vcf` / `lab_b.vcf`，
以及一份含非法 REF 的 `lab_bad.vcf` 用于演示失败）。

## 比对规则

**记录合法性**（任一记录非法 → 整次比对失败，失败原因含文件与行号）：

- 每份 VCF 数据记录 ≤ 1000 条；染色体名须与 FASTA 序列名一致；
- POS 为正整数且不越出参考序列；REF 与参考序列逐碱基核对一致；
- 仅接受 A/C/G/T 的双等位 SNV 或纯插入缺失：REF≠ALT，ALT 不得为多等位（`,`）、
  符号型（`<DEL>`）、`*` 等；REF/ALT 含 N 或其他字符亦非法；
- 去除共有前后缀后两端都非空的记录（MNP、复杂替换）视为非法。

**归一化（等价写法归为同一位点）**：

1. 核对 REF 后，先去除 REF/ALT 的共有后缀、再去共有前缀；
2. 纯插入 / 纯缺失按参考序列**左对齐**（在串联重复、均聚物中滑到最左代表）；
3. 归一化后的 `(pos, ref, alt)` 为位点唯一键：插入的 `ref` 为空（表示在 pos 之前插入），
   缺失的 `alt` 为空，网页中以 `-` 显示。

**合并与归类**：同一文件内归一化键相同的记录合并为一个位点，保留全部原行号
（行号按上传文件的物理行计，含表头）；两份文件按键求交/差，分为共有、甲独有、乙独有。

## 可靠性与并发

- 每次上传在一个事务内写入自己的报告与全部位点（`reports` / `variants` 表，
  按 `report_id` 关联），失败整体回滚，不留半成品；
- 非法记录只产生一条"失败"报告，**不改动任何旧报告**；
- 处理过程全部使用请求内的局部变量与每请求独立的数据库连接（WAL 模式），
  并发上传的报告不会混用不同批次的参考序列、位点和行号。

## 目录结构

```
app/
  compare.py    # FASTA/VCF 解析、REF 核对、左对齐归一化、集合比较（纯逻辑）
  db.py         # SQLite 存取（WAL、单事务写入）
  main.py       # FastAPI 路由与页面渲染
  templates/    # 上传表单、报告列表、报告详情
tests/          # 单元测试 + 网页集成测试 + 并发隔离测试
examples/       # 示例参考序列与 VCF
Dockerfile / docker-compose.yml
```
