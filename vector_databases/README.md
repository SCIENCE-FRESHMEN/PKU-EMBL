# CAZySeek 向量数据库 · 总说明（README）

> 本目录是 **CAZySeek** 项目的交付核心：四个独立、可检索的 **ChromaDB 向量数据库**，作为后续构建**层级式 GraphRAG（Hierarchical GraphRAG）** 的数据底座。
> 面向：北大环境微生物生物技术课题组 · CAZyme（碳水化合物活性酶）生物挖掘（Bio-mining）。
> 工程标准对齐 *ReactionSeek*（**Nature Communications** 2026, doi:10.1038/s41467-026-70180-1）。

---

## 0. 一句话概览

- **4 个库 = 4 类实体**：概念（CAZypedia）/ 酶家族（CAZy）/ 反应边（Substrate）/ 基因簇（CGC）。
- **共享 1 个嵌入模型**：`all-MiniLM-L6-v2`（384 维，cosine）。
- **双组件架构**（对齐论文）：每个库 = 向量组件（ChromaDB）+ 补充元数据组件（`supplementary.json`）。
- **层级主干键**：`family_id`（01/02/04）与 `cazyme_families_base`（03 的列表型字段），把四库串成 `Class → Family → Reaction → Cluster`。

---

## 一、目录结构与「每个文件夹是什么 / 怎么用」

```
vector_databases/
├── 01_cazypedia/             # 库①：CAZypedia 百科（概念节点）
├── 02_cazy_database/         # 库②：CAZy 家族统计与活性（酶家族节点）
├── 03_cgc_database/          # 库③：dbCAN-seq 基因簇 + 预测底物（基因簇节点）
├── 04_substrate_specificity/ # 库④：CAZy 活性反应表（反应关系边）
├── models/
│   └── all-MiniLM-L6-v2/     # 四库共用的嵌入模型（88 MB）
└── README.md                 # 本文件
```

### 1.1 每个文件夹的用途与打开方式

| 文件夹 | 是什么 | 体量 | 怎么用（打开） |
|--------|--------|------|----------------|
| `01_cazypedia/` | CAZypedia 叙述性概念/机制文本分块 | 1,877 条 · 24 MB | `chromadb.PersistentClient(path="01_cazypedia").get_collection("cazypedia")` |
| `02_cazy_database/` | cazy.org 六类家族（GH/GT/PL/CE/AA/CBM）统计与已表征活性 | 1,618 条 · 12 MB | `...get_collection("cazy_db")` |
| `03_cgc_database/` | dbCAN-seq 预测的 CAZyme 基因簇（CGC），含家族组成、预测底物、物种 | 12,121 条 · 79 MB | `...get_collection("cgc_db")` |
| `04_substrate_specificity/` | 来自 CAZy 活性表的结构化反应记录（底物/产物/酶/EC） | 862 条 · 5.6 MB | `...get_collection("substrate_specificity")` |
| `models/all-MiniLM-L6-v2/` | 共享嵌入模型 | 88 MB | 由 Sentence-Transformers 直接加载，路径见下 |

> **合计：15,578 条向量，约 205 MB（含模型）。**

### 1.2 每个库目录的内部结构（四个库一致）

```
<库目录>/
├── chroma.sqlite3      # ChromaDB 主库：集合元数据 + 向量索引指针
├── <uuid>/             # HNSW 向量索引段（实际 384 维向量存放处）
│   ├── index_<n>.bin
│   └── data_<n>.bin
└── supplementary.json  # 补充元数据仓库（双组件架构的「第二组件」）
```

**双组件怎么用：**
- **向量组件（chroma.sqlite3 + HNSW 段）**：负责相似度检索，由 ChromaDB 管理，无需手动解析。
- **补充元数据组件（supplementary.json）**：以 `chunk_id` 为键，保存该条**完整文本 + 全部结构化字段**（家族号、EC、底物、物种等）。检索命中后，用返回的 `id` 在 `supplementary.json` 中取回原文与细节——避免把长文本全塞进向量库元数据。
- 二者通过 `chunk_id` / `id` 一一对应，数量始终一致（`verify_all.py` 会校验）。

### 1.3 最小可用代码（加载模型 + 打开任意一个库）

```python
import chromadb, os, json
from sentence_transformers import SentenceTransformer

HERE = os.path.dirname(os.path.abspath(__file__))          # 本 vector_databases/ 目录
MODEL = SentenceTransformer(os.path.join(HERE, "models", "all-MiniLM-L6-v2"))

def open_db(subdir, coll):
    client = chromadb.PersistentClient(path=os.path.join(HERE, subdir))
    return client.get_collection(coll)

def supp(subdir):
    return json.load(open(os.path.join(HERE, subdir, "supplementary.json"), encoding="utf-8"))
```

---

## 二、每个数据库里的「表」是什么 / 怎么用

> ChromaDB 中每个库只有 **一个集合（collection）**，语义上等同于一张「表」。
> 该表的「列」由四部分构成：**`ids`**（行主键）、**`documents`**（可检索文本）、**`embeddings`**（384 维向量，通常不手动读）、**`metadatas`**（结构化过滤/路由字段）。
> 下面按库给出 `metadatas` 的列定义、含义与典型用法。

### 表 ① `cazypedia`（库 01_cazypedia）

| 列（metadata 字段） | 类型 | 含义 | 用法 |
|---------------------|------|------|------|
| `class` | str | CAZy 大类（GH/GT/PL/CE/AA/CBM/Lexicon） | 按类路由：`where={"class": "$eq": "GH"}` |
| `family_id` | str | 家族号，如 `AA10`（概念/词表页可能为空） | 跨库关联键，与 02/04 同号 |
| `level` | str | 层级：`family` / `concept` 等 | 区分「家族页」与「概念页」 |
| `section` | str | 维基分节标题（如 `Introduction`） | 定位文本出处 |
| `title` | str | 页面标题 | 展示用 |
| `url` | str | 来源 URL | 溯源 |

- **`documents`**：CAZypedia 词条按 `== 章节 ==` 切分后的文本块。
- **典型查询**：
  ```python
  col = open_db("01_cazypedia", "cazypedia")
  col.query(query_texts=["cellulase hydrolysis of cellulose"], n_results=5)          # 语义检索
  col.get(where={"family_id": {"$eq": "GH1"}})                                       # 精确家族路由
  ```

### 表 ② `cazy_db`（库 02_cazy_database）

| 列（metadata 字段） | 类型 | 含义 | 用法 |
|---------------------|------|------|------|
| `class` | str | CAZy 大类 | 类级路由 |
| `family_id` | str | 家族号，如 `AA0`、`GH1` | **核心关联键**（`$eq` 精确匹配） |
| `section` | str | 固定为 `family_summary` | 标识文本类型 |
| `source` | str | 固定为 `CAZy_Database` | 来源标记 |
| `url` | str | 家族页 URL（如 `https://www.cazy.org/GH1.html`） | 溯源 |

- **`documents`**：家族页摘要文本（家族定义、统计数、已表征活性概述）。
- **典型查询**：
  ```python
  col = open_db("02_cazy_database", "cazy_db")
  col.query(query_texts=["GH1 beta-glucosidase characterized activity"], n_results=3)
  col.get(where={"family_id": {"$eq": "GH1"}})     # 精确命中酶家族节点
  ```

### 表 ③ `cgc_db`（库 03_cgc_database）

| 列（metadata 字段） | 类型 | 含义 | 用法 |
|---------------------|------|------|------|
| `cgc_id` | str | 基因簇 ID，如 `MGYG000290072_56\|CGC1` | 行主键补充 |
| `genome_id` | str | 基因组 ID（MAG/基因组的唯一标识） | 回查物种/GTDB 补注 |
| `substrate` | str | dbCAN-seq 预测的该簇底物（如 `agarose`、`xylan`） | 按底物检索基因簇 |
| `species` | str | 物种（约 22% 空缺，源自 MAG 未注释） | 可选过滤 |
| `cazyme_families_base` | **list** | 该簇所含的 CAZyme **基础家族**组合，如 `["GH16","GH43","CBM32"]` | **层级关联键**，`$contains "GH1"` 全词匹配 |
| `cazyme_families` | **list** | 含亚家族后缀的家族，如 `["GH43_19","GH43_24"]` | 更细粒度家族定位 |
| `classes` | **list** | 由家族推导出的类列表，如 `["CBM","GH"]` | 类级路由 |
| `num_cazyme` | int | 簇内 CAZyme 数量 | 规模过滤 |
| `num_genes` | int | 簇内基因总数 | 规模过滤 |
| `source` | str | 固定为 `CGC_Database` | 来源标记 |

- **`documents`**：该 CGC 的自然语言描述（基因组 + 家族组成 + 预测底物 + 物种）。
- **典型查询**：
  ```python
  col = open_db("03_cgc_database", "cgc_db")
  # 语义：找作用于木聚糖的基因簇
  col.query(query_texts=["gene cluster acting on xylan substrate"], n_results=5)
  # 精确（全词）：找所有含 GH1 的基因簇 —— 不会误命中 GH13/GH16
  col.get(where={"cazyme_families_base": {"$contains": "GH1"}})
  # 类级路由
  col.get(where={"classes": {"$contains": "GH"}})
  ```

### 表 ④ `substrate_specificity`（库 04_substrate_specificity）

| 列（metadata 字段） | 类型 | 含义 | 用法 |
|---------------------|------|------|------|
| `activity_name` | str | 反应名称，如 `lytic chitin monooxygenase` | 文本检索 |
| `class` | str | CAZy 大类 | 类级路由 |
| `ec` | str | EC 编号，如 `1.14.99.53` | 按酶学编号定位 |
| `family_id` | str | 家族号，如 `AA10` | **核心关联键**（`$eq`） |
| `source` | str | 固定为 `Substrate_Specificity` | 来源标记 |
| `url` | str | 家族页 URL | 溯源 |

- **`documents`**：结构化反应描述（底物 / 产物 / 作用键 / EC / 反应类型）。
- **典型查询**：
  ```python
  col = open_db("04_substrate_specificity", "substrate_specificity")
  col.query(query_texts=["sialic acid cleavage by neuraminidase"], n_results=3)
  col.get(where={"family_id": {"$eq": "GH1"}})     # 取该家族的全部反应边
  col.get(where={"ec": {"$eq": "3.2.1.23"}})       # 按 EC 号精确定位
  ```

---

## 三、跨库联合检索（代码示例）

```python
import chromadb, os
HERE = os.path.dirname(os.path.abspath(__file__))
DBS = {                                   # collection -> 目录
    "cazypedia":            "01_cazypedia",
    "cazy_db":              "02_cazy_database",
    "cgc_db":               "03_cgc_database",
    "substrate_specificity":"04_substrate_specificity",
}
_COLS = {c: chromadb.PersistentClient(path=os.path.join(HERE, p)).get_collection(c)
         for c, p in DBS.items()}

def search_all(query, top_k=3):
    """四库语义联合检索，返回 {collection: (documents, metadatas, ids)}"""
    out = {}
    for coll, col in _COLS.items():
        r = col.query(query_texts=[query], n_results=top_k)
        out[coll] = (r["documents"][0], r["metadatas"][0], r["ids"][0])
    return out

# 例：一次问询，四库分别命中
res = search_all("β-glucosidase that hydrolyzes cellobiose", top_k=2)
```

> 注：若查询时不显式传嵌入，ChromaDB 会调用其内置的 all-MiniLM-L6-v2 ONNX 实现，与本地 Sentence-Transformers 输出等价，结果一致。

---

## 四、面向「层级 RAG（Hierarchical GraphRAG）」怎么用本数据库

### 4.1 四库 → 四类图节点 + 关系边

| 图元素 | 来源库 | 节点/边 | 定位方式 |
|--------|--------|---------|----------|
| **概念节点** | 01 CAZypedia | 节点 | 语义检索 + `family_id`/`class` 路由 |
| **酶家族节点** | 02 CAZy | 节点 | `family_id $eq`（精确） |
| **反应关系边** | 04 Substrate | 边（家族→底物/产物） | `family_id $eq` |
| **基因簇节点** | 03 CGC | 节点 | `cazyme_families_base $contains`（全词） |

### 4.2 层级主干（共享键串接）

```
        Class (GH/GT/PL/CE/AA/CBM)
          │  family_id / class
          ▼
        Family (GH1)  ── 02 CAZy 节点
          │  ├─ family_id $eq ── 04 Substrate 反应边
          │  └─ family_id $eq ── 01 CAZypedia 概念注释
          │  cazyme_families_base $contains ── 03 CGC 基因簇节点
          ▼
   Subfamily / Cluster (GH43_19, MGYG…_CGC1)
```

- **01 ↔ 02 ↔ 04** 通过标量 `family_id`（`$eq`）直接连边；
- **02 ↔ 03** 通过 03 的**列表型** `cazyme_families_base`（`$contains` 全词）连边；
- 由此可自动建边：**概念（01）— 酶家族（02）— 反应（04）— 基因簇（03）**，形成层级主干。

### 4.3 推荐的检索 / 路由流程（供 GraphRAG 检索器复用）

1. **类路由（Class Router）**：用 `class` 字段把查询分派到相关大类（GH/GT/…）。
2. **家族路由（Family Router）**：抽取查询中的家族号（如 `GH1`），用 `$eq` 在 01/02/04 精确命中，用 `$contains` 在 03 命中含该家族的基因簇。
3. **语义兜底**：家族号缺失时，用 `query_texts` 做四库语义检索（`search_all`）。
4. **实体抽取输入**：检索命中后，用 `id` 到对应 `supplementary.json` 取回**完整结构化字段**，直接作为后续 LLM 实体抽取 / 关系建边的上下文。

### 4.4 实测贯通示例（以 `GH1` 为主干，ChromaDB 检索层）

| 关联 | 查询 | 命中 |
|------|------|------|
| 02 酶家族节点 | `family_id $eq GH1` | 1 条 |
| 04 反应边 | `family_id $eq GH1` | 3 条（xylan exo-xylosidase、arabinan exo-arabinofuranosidase 等） |
| 03 基因簇节点 | `cazyme_families_base $contains GH1` | 5 条，且均真实含全词 `GH1`（无误命中 GH13/GH16） |

→ **Class → Family → Reaction → Cluster 主干在检索层可直接贯通**，可作为层级 GraphRAG 的实体定位与建边基础设施。

---

## 五、与参考论文（ReactionSeek）标准完成对比

对照基线：论文 **Design of SynChat** 章节的明文工程标准（已从 PDF 抽取原文核对）。

| 编号 | 论文标准 | 论文原文（摘） | 本库完成情况 |
|------|----------|----------------|--------------|
| **S1** | 向量库 = ChromaDB | *"indexed and managed using **ChromaDB**, an open-source vector database"* | ✅ 四库均为 ChromaDB 持久化实例 |
| **S2** | 嵌入 = all-MiniLM-L6-v2（SBERT，384 维，无降维） | *"**all-MiniLM-L6-v2** … **384-dimensional dense vector embeddings without dimensionality reduction**"* | ✅ 四库共用该模型；维度复核 = 384，归一化稠密向量 |
| **S3** | 双组件：向量查询库 + 补充元数据仓库（JSON） | *"a vectorized query database … and a **supplementary metadata repository**"* | ✅ 每库均含 `chroma.sqlite3` + `supplementary.json`，条数一致 |
| **S4** | 检索 = 向量语义相似度 + ChromaDB 内**精确字符匹配（全词）** | *"**vector similarity (semantic match) and exact character matching within the ChromaDB index**"* | ✅ 语义检索 + `$eq`/`$contains` 精确匹配；03 全词反例已验证（GH1 不误命中 GH13/GH16） |
| **S5** | 距离度量 = Cosine | Fig.8 *"Cosine similarity"* | ✅ 建库 `hnsw:space = "cosine"` |
| **S6** | 层级 metadata（支撑层级 RAG 路由） | 论文上溯「反应文档→参数」层级 | ✅ 各库带层级键：01/02/04 用 `family_id`/`class`，03 用 `cazyme_families_base`(list)/`classes`(list) |

> **结论：四库对论文 S1–S6 六项明文标准全部达标（✅）。**

---

## 六、已知限制与后续建议（非阻断）

| 项 | 说明 | 建议 |
|----|------|------|
| **04 仅覆盖 231/525 家族** | 数据本质：CAZy 只对有「已表征底物特异性」的家族维护结构化反应表（以 GH/CE/PL 为主）；CBM 非催化、AA 氧化还原、未策展 GT/CE/PL 自然缺失 | 若需补 GT 等反应边，后续接入 CAZy `Activities` 或 UniProt/文献反应注释 |
| **03 物种字段 22% 空缺** | 源自 dbCAN-seq MAG 未做物种级注释（仅 genome_id） | 用 `genome_id` 回查 GTDB/NCBI taxonomy 补注「基因簇—物种」边 |
| **01 概念页 `family_id` 可能为空** | 词表/概念页不绑定单一家族 | GraphRAG 建边时以非空 `family_id` 页为主干，概念页作注释层 |
| **03 体量最大且持续增长** | 12,121 条，随 dbCAN-seq 扩充 | 建议独立生命周期管理（月度/按需重建），其余库季度更新 |

---

*构建标准对齐 ReactionSeek（NC 2026）。嵌入模型 all-MiniLM-L6-v2，cosine 距离，ChromaDB 持久化。*
