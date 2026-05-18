[README.md](https://github.com/user-attachments/files/27951063/README.md)
# Protein Linker Analysis Toolkit

蛋白质Linker结构分析工具集，用于分析RFdiffusion生成的PDB文件中linker结构的几何参数，评估蛋白质刚柔性，筛选优质序列。

## 功能特性

### 1. Linker结构分析 (`analyze_linker_pdbs.py`)
- **PDB解析**：自动解析PDB文件中的原子坐标信息
- **Linker检测**：识别A链上的连续GLY残基作为linker区域
- **几何分析**：
  - CA-CA步长偏差计算
  - 肽键(C-N)长度检测
  - 二面角(ω)分析
  - 回转半径(Rg)计算
  - 端到端距离比(L/L0)
- **冲突检测**：
  - 链内空间冲突检测
  - A-B链间冲突检测
- **评分排序**：基于多维度指标的综合评分系统
- **输出**：CSV和Word格式的分析报告

### 2. Linker序列比对 (`linker_alignment.py`)
- 从完整序列中提取linker区域
- 多序列比对可视化
- 支持对照序列(cIgg2)比较
- 生成PNG格式的比对图和FASTA序列文件
- 氨基酸颜色编码（按性质分类）

### 3. 蛋白质综合评分 (`protein_score.py`)
- **表达性评分**：基于GC含量、稀有氨基酸、序列复杂度
- **水溶性评分**：亲水/疏水比例、净电荷、等电点
- **结构评分**：脯氨酸、甘氨酸比例、二级结构倾向
- **稳定性评分**：二硫键潜力、疏水核心、稳定性氨基酸
- **综合评分**：加权平均（表达性20%、水溶性25%、结构25%、稳定性30%）

### 4. 刚柔性分析 (`flexibility_analysis.py`)
- 柔性/刚性氨基酸比例分析
- Pro-Gly motif检测
- 疏水核心评估
- 二硫键潜力分析
- 综合评分：表达性、柔性、稳定性、平衡性

### 5. A链序列提取 (`python extract_A_chain.py`)
- 从PDB文件中批量提取A链氨基酸序列
- 支持FASTA格式输出
- 自动去重重复残基

## 安装

```bash
pip install python-docx scipy numpy matplotlib
```

## 使用方法

### 1. Linker结构分析
```bash
python analyze_linker_pdbs.py
# 输出:
# - linker_analysis_ranking_strict_with_B.csv
# - linker_analysis_report_strict_zh_with_B.docx
```

### 2. 序列比对
```bash
python linker_alignment.py
# 需要提前准备 top50_sequences.fasta 文件
# 输出:
# - linker_alignment.png (可视化比对图)
# - linker_alignment_sequences.fasta
```

### 3. 蛋白质评分
```bash
python protein_score.py
# 需要提前准备 chain_a_sequences.fasta 文件
# 输出:
# - protein_scores.txt (评分表)
# - top50_sequences.fasta (前50名序列)
```

### 4. 刚柔性分析
```bash
python flexibility_analysis.py
# 需要提前准备 protein_scores.txt 和 chain_a_sequences.fasta
# 输出:
# - flexibility_analysis.txt
# - top50_flexible_sequences.fasta
```

### 5. 提取A链序列
```bash
python "python extract_A_chain.py"
# 输出: A链序列合集.txt
```

## 项目结构

```
analyze_linker_tool/
├── analyze_linker_pdbs.py       # Linker结构分析
├── linker_alignment.py          # 序列比对可视化
├── protein_score.py            # 蛋白质综合评分
├── flexibility_analysis.py      # 刚柔性分析
├── python extract_A_chain.py    # A链序列提取
├── README.md                   # 本文件
├── LICENSE                     # MIT许可证
└── requirements.txt           # 依赖包列表
```

## 依赖

- Python 3.8+
- python-docx >= 0.8.10
- scipy >= 1.7.0
- numpy >= 1.21.0
- matplotlib >= 3.4.0 (用于序列比对可视化)

## 评分标准

### Linker分析评分

| 指标 | 说明 | 扣分范围 |
|------|------|----------|
| CA-CA步长 | 理想值3.8A，合理区间2.8-4.5A | 0-38分 |
| 肽键长度 | 理想值1.33A，<1.15压缩，>2.2断裂 | 0-60分 |
| 二面角ω | trans≈180°，cis≈0° | 0-20分 |
| 空间冲突 | <2.0A严重，<2.4A轻度 | 0-50分 |
| 形状指标 | Rg和L/L0(柔性linker常见0.4-0.6) | 0-24分 |

### 刚柔性分析理想参数

| 参数 | 理想范围 |
|------|----------|
| 柔性氨基酸比例 | 25-40% |
| 刚性氨基酸比例 | 20-35% |
| 疏水核心比例 | 20-35% |
| 柔刚比 | 0.7-1.5 |
| 甘氨酸比例 | 3-10% |
| 脯氨酸比例 | 3-8% |
| 二硫键数量 | >=2 |

## 许可证

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

## 作者

[Your Name]

## 贡献

欢迎提交Issue和Pull Request！
