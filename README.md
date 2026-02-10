# plantfamilyallin

**plantfamilyallin** 是一个基于 Snakemake 的植物基因家族系统分析流程，
用于从全基因组水平系统解析植物基因家族的鉴定、结构特征、系统发育、
共线性关系及选择压力（Ka/Ks）。

该流程支持多物种比较，具备良好的可重复性与可扩展性，
适用于植物功能基因组学与进化基因组学研究。

---

## ✨ 功能模块概览

1. **基因组注释清洗与最长转录本提取**
2. **基因家族成员鉴定**
   - BLAST
   - Pfam 结构域
   - HMM 搜索
3. **染色体定位与分布可视化**
4. **蛋白 Motif（MEME）与基因结构分析**
5. **启动子序列提取与顺式作用元件分析**
6. **蛋白理化性质与亚细胞定位预测**
7. **多物种共线性分析（MCScanX / GENESPACE）**
8. **系统发育树构建与可视化**
9. **家族内 Ka/Ks 选择压力分析**
10. **全基因组共线锚点 Ks 背景分布对比**

---

## 🧩 软件架构

```text
plantfamilyallin/
├── snakfile                 # 主流程 Snakefile
├── config.yaml              # 用户配置文件
├── envs/
│   └── plantfamilyallin.yaml
├── scripts/                 # Python / R / Shell 工具脚本
├── tools/
│   └── KaKs_Calculator-3.0/ # 外部工具（需用户自行准备）
└── results/                 # 输出结果目录


软件环境要求
操作系统

Linux（推荐 Ubuntu / CentOS）

硬件建议

CPU：≥ 8 cores

内存：≥ 32 GB（大基因组建议 ≥ 64 GB）

磁盘：≥ 100 GB（视物种数量而定）


依赖环境安装（Conda）
1安装 Conda（如未安装）
https://docs.conda.io/en/latest/miniconda.html

2️创建运行环境
conda env create -f envs/plantfamilyallin.yaml
conda activate plantfamilyallin

3 检查关键软件
snakemake --version
diamond --version
hmmscan -h
mcscanx -h
iqtree2 -h

🔧 外部依赖：KaKs_Calculator

Ka/Ks 分析模块使用 KaKs_Calculator-3.0，
该工具不通过 Conda 安装，需用户自行准备。

1️ 下载
http://ngdc.cncb.ac.cn/biocode/tools/BT000001

2️ 放置目录结构
tools/
└── KaKs_Calculator-3.0/
    ├── bin/
    │   ├── KaKs
    │   └── AXTConvertor
    └── pal2nal.pl

3️ 在 config.yaml 中指定路径
kaks:
  kaks_bin_dir: "tools/KaKs_Calculator-3.0/bin"
  pal2nal: "tools/KaKs_Calculator-3.0/pal2nal.pl"

📝 配置文件说明（config.yaml）

核心参数示例：

family_name: "bHLH"

target:
  name: "SL"
  genome_fa: "data/SL/SL.fasta"
  gff3: "data/SL/SL.gff3"

model_family_pep: "data/model/AT_bHLH.fa"

pfam_hmm: "data/Pfam/PF00010.hmm"
pfam_domains_of_interest:
  - PF00010

threads: 10
outdir: "results"


支持多物种共线分析：

synteny_species:
  - name: "AT"
    genome_fa: "data/AT/TAIR10.fa"
    gff3: "data/AT/TAIR10.gff3"

▶️ 运行流程

在项目根目录执行：

snakemake -j 10 --use-conda


或后台运行：

nohup snakemake -j 20 --use-conda > run.log 2>&1 &

📂 输出结果说明
results/
├── 01.cds_protein/
├── 02.family_id/
├── 03.chromosome_map/
├── 04.meme_structure/
├── 05.promoter_cis/
├── 06.protein_property/
├── 07.synteny/
├── 08.phylogeny/
├── 09.selection/
├── 10.syntenic_kaks/
└── 99.result/   # 所有最终 PDF 图件

📊 主要输出图件

家族鉴定 Venn 图

染色体定位图

Motif + Gene Structure 综合图

系统发育树

Ka/Ks 分布图

家族 vs 全基因组 Ks 对比图

🔁 可重复性说明

使用 Snakemake 管理工作流

提供 Conda 环境定义文件

所有路径均为相对路径

支持容器化部署（Docker / Apptainer）

📄 许可证

本软件仅用于科研用途，
如用于商业用途请联系作者。

📮 联系方式

作者：JWJ
邮箱：13164328557@163.com





