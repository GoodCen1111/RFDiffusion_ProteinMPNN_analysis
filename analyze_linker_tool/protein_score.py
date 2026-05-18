"""
蛋白质可表达性、水溶性、结构和稳定性综合评估脚本
"""

import os
import re
from collections import OrderedDict

# 氨基酸理化性质
AA_PROPERTIES = {
    'A': {'hydrophobicity': 1.8, 'charge': 0, 'volume': 88.6, 'polar': False},
    'R': {'hydrophobicity': -4.5, 'charge': 1, 'volume': 173.4, 'polar': True},
    'N': {'hydrophobicity': -3.5, 'charge': 0, 'volume': 114.1, 'polar': True},
    'D': {'hydrophobicity': -3.5, 'charge': -1, 'volume': 111.1, 'polar': True},
    'C': {'hydrophobicity': 2.5, 'charge': 0, 'volume': 108.5, 'polar': False},
    'E': {'hydrophobicity': -3.5, 'charge': -1, 'volume': 138.4, 'polar': True},
    'Q': {'hydrophobicity': -3.5, 'charge': 0, 'volume': 143.8, 'polar': True},
    'G': {'hydrophobicity': -0.4, 'charge': 0, 'volume': 60.1, 'polar': False},
    'H': {'hydrophobicity': -3.2, 'charge': 0.5, 'volume': 153.2, 'polar': True},
    'I': {'hydrophobicity': 4.5, 'charge': 0, 'volume': 166.7, 'polar': False},
    'L': {'hydrophobicity': 3.8, 'charge': 0, 'volume': 166.7, 'polar': False},
    'K': {'hydrophobicity': -3.9, 'charge': 1, 'volume': 168.6, 'polar': True},
    'M': {'hydrophobicity': 1.9, 'charge': 0, 'volume': 162.9, 'polar': False},
    'F': {'hydrophobicity': 2.8, 'charge': 0, 'volume': 189.9, 'polar': False},
    'P': {'hydrophobicity': -1.6, 'charge': 0, 'volume': 112.7, 'polar': False},
    'S': {'hydrophobicity': -0.8, 'charge': 0, 'volume': 89.0, 'polar': True},
    'T': {'hydrophobicity': -0.7, 'charge': 0, 'volume': 116.1, 'polar': True},
    'W': {'hydrophobicity': -0.9, 'charge': 0, 'volume': 227.8, 'polar': False},
    'Y': {'hydrophobicity': -1.3, 'charge': 0, 'volume': 193.6, 'polar': False},
    'V': {'hydrophobicity': 4.2, 'charge': 0, 'volume': 140.0, 'polar': False},
}

# 稀有氨基酸（在大肠杆菌中低频使用）
RARE_AA = set('RWQMEKFYI')
# 常用氨基酸
COMMON_AA = set('AGSTNDEPK')

# 脯氨酸和甘氨酸（对折叠有特殊影响）
PRO_GLY = set('PG')
# 半胱氨酸（形成二硫键）
CYS = set('C')
# 疏水氨基酸
HYDROPHOBIC = set('AILFWVMP')
# 亲水氨基酸
HYDROPHILIC = set('RKDENQSTYH')


def read_fasta(fasta_file):
    """读取FASTA文件"""
    sequences = OrderedDict()
    current_id = None
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)

        if current_id:
            sequences[current_id] = ''.join(current_seq)

    return sequences


def calculate_composition(seq):
    """计算氨基酸组成"""
    length = len(seq)
    comp = {}
    for aa in AA_PROPERTIES.keys():
        count = seq.count(aa)
        comp[aa] = count / length * 100 if length > 0 else 0
    return comp


def calculate_expressibility_score(seq, comp):
    """
    表达性评分 (0-100)
    - GC含量模拟（基于疏水性氨基酸）
    - 稀有氨基酸比例
    - 序列复杂度
    """
    length = len(seq)
    if length == 0:
        return 0

    score = 100

    # 稀有氨基酸惩罚（稀有氨基酸过多不利于表达）
    rare_count = sum(comp[aa] for aa in RARE_AA)
    score -= rare_count * 2

    # 常用氨基酸加分
    common_count = sum(comp[aa] for aa in COMMON_AA)
    score += common_count * 0.5

    # 序列长度惩罚（过长或过短都不利于表达）
    if length < 50:
        score -= (50 - length) * 0.5
    elif length > 500:
        score -= (length - 500) * 0.1

    # 连续的相同氨基酸惩罚
    repeats = max([len(m.group(0)) for m in re.finditer(r'(.)\1{2,}', seq)] or [0])
    score -= repeats * 2

    # 天冬酰胺(N)和谷氨酰胺(Q)比例（容易发生脱酰胺）
    nq_ratio = (comp['N'] + comp['Q']) / 100
    if nq_ratio > 0.15:
        score -= (nq_ratio - 0.15) * 50

    return max(0, min(100, score))


def calculate_solubility_score(seq, comp):
    """
    水溶性评分 (0-100)
    - 亲水/疏水比例
    - 净电荷
    - 等电点预测
    """
    length = len(seq)
    if length == 0:
        return 0

    score = 100

    # 亲水性评估（基于氨基酸的亲疏水性）
    hydro_score = sum(AA_PROPERTIES[aa]['hydrophobicity'] for aa in seq if aa in AA_PROPERTIES)
    avg_hydro = hydro_score / length

    # 理想蛋白质的平均疏水性约在-1到1之间
    if avg_hydro > 1.5:
        score -= (avg_hydro - 1.5) * 20  # 太疏水
    elif avg_hydro < -2:
        score -= (-2 - avg_hydro) * 10  # 太亲水

    # 疏水氨基酸比例（理想范围20-35%）
    hydro_ratio = sum(comp[aa] for aa in HYDROPHOBIC) / 100
    if hydro_ratio < 0.20:
        score -= (0.20 - hydro_ratio) * 50
    elif hydro_ratio > 0.40:
        score -= (hydro_ratio - 0.40) * 30
    elif 0.20 <= hydro_ratio <= 0.35:
        score += 10  # 理想疏水比例加分

    # 亲水氨基酸加分
    hydroPhilic_ratio = sum(comp[aa] for aa in HYDROPHILIC) / 100
    if hydroPhilic_ratio > 0.30:
        score += (hydroPhilic_ratio - 0.30) * 20

    # 净电荷（带电氨基酸多有利于水溶性）
    net_charge = sum(AA_PROPERTIES[aa]['charge'] for aa in seq if aa in AA_PROPERTIES)
    charge_score = abs(net_charge) / length * 100
    score += min(charge_score, 15)

    # 碱性氨基酸（K,R,H）有利于水溶性
    basic_ratio = (comp['K'] + comp['R'] + comp['H']) / 100
    if basic_ratio > 0.08:
        score += (basic_ratio - 0.08) * 30

    # 酸性氨基酸（E,D）也有利于水溶性
    acidic_ratio = (comp['E'] + comp['D']) / 100
    if acidic_ratio > 0.08:
        score += (acidic_ratio - 0.08) * 20

    return max(0, min(100, score))


def calculate_structure_score(seq, comp):
    """
    结构合理性评分 (0-100)
    - 脯氨酸和甘氨酸比例
    - 二级结构倾向
    - 无序区域检测
    """
    length = len(seq)
    if length == 0:
        return 0

    score = 100

    # 脯氨酸比例（脯氨酸多会破坏二级结构，但适量有助于折叠）
    pro_ratio = comp['P'] / 100
    if pro_ratio > 0.08:
        score -= (pro_ratio - 0.08) * 200
    elif pro_ratio < 0.02:
        score -= 5  # 脯氨酸太少也不利于结构

    # 甘氨酸比例（甘氨酸多增加柔性，但过多导致无序）
    gly_ratio = comp['G'] / 100
    if gly_ratio > 0.10:
        score -= (gly_ratio - 0.10) * 150
    elif gly_ratio < 0.03:
        score -= 5

    # Pro-Gly序列（常见于转角）
    pg_count = len(re.findall(r'P.G|G.P', seq))
    score += min(pg_count * 2, 10)

    # 半胱氨酸比例（可形成二硫键）
    cys_ratio = comp['C'] / 100
    if cys_ratio > 0.02:
        # 理想的半胱氨酸比例（可形成适当的二硫键）
        score += min(cys_ratio * 200, 15)

    # 疏水核心评估（应该有交替的疏水残基）
    hydrophobic_stretches = len(re.findall(r'[AILFWVMP]{5,}', seq))
    if hydrophobic_stretches > 0:
        score -= hydrophobic_stretches * 5  # 连续疏水太长不利于折叠

    # 芳香族氨基酸比例（影响蛋白质稳定性）
    aromatic = comp['F'] + comp['Y'] + comp['W']
    if aromatic > 0.08:
        score += (aromatic - 0.08) * 30  # 适当芳香族增加稳定性

    # 正负电荷残基的分布（电荷促进折叠）
    charged_pattern = len(re.findall(r'[RK][DE]|[DE][RK]', seq))
    score += min(charged_pattern, 10)

    # 丝氨酸和苏氨酸比例（磷酸化位点，但也有利于水溶性）
    st_ratio = (comp['S'] + comp['T']) / 100
    if 0.05 <= st_ratio <= 0.12:
        score += 5

    return max(0, min(100, score))


def calculate_stability_score(seq, comp):
    """
    稳定性评分 (0-100)
    - 二硫键潜力
    - 疏水核心
    - 稳定性氨基酸
    """
    length = len(seq)
    if length == 0:
        return 0

    score = 100

    # 半胱氨酸数量（形成二硫键的关键）
    cys_count = seq.count('C')
    if cys_count >= 2:
        # 估算可能形成的二硫键数量
        potential_disulfide = cys_count // 2
        score += min(potential_disulfide * 5, 20)
    elif cys_count == 1:
        score -= 5  # 只有一个半胱氨酸不能形成二硫键

    # 疏水核心稳定性（疏水氨基酸应分布在整个序列中）
    hydrophobic_count = sum(1 for aa in seq if aa in HYDROPHOBIC)
    hydro_ratio = hydrophobic_count / length

    # 理想的疏水比例（约25-35%）
    if 0.25 <= hydro_ratio <= 0.35:
        score += 10
    elif hydro_ratio < 0.20:
        score -= (0.20 - hydro_ratio) * 30
    elif hydro_ratio > 0.40:
        score -= (hydro_ratio - 0.40) * 20

    # 稳定性氨基酸加分（缬氨酸、异亮氨酸亮氨酸等）
    stabilizing = comp['V'] + comp['I'] + comp['L'] + comp['M']
    score += min(stabilizing * 0.3, 10)

    # 不稳定性氨基酸惩罚（天冬酰胺、谷氨酰胺容易脱酰胺）
    unstable_aa = comp['N'] + comp['Q']
    if unstable_aa > 15:
        score -= (unstable_aa - 15) * 2

    # 脯氨酸的顺反异构（脯氨酸多可能影响稳定性）
    if comp['P'] > 8:
        score -= (comp['P'] - 8) * 3

    # 甘氨酸过多导致柔性过高（不稳定）
    if comp['G'] > 10:
        score -= (comp['G'] - 10) * 3

    # 甲硫氨酸氧化（Met在表面容易被氧化）
    met_ratio = comp['M'] / 100
    if met_ratio > 0.02:
        score -= met_ratio * 50

    # 带电氨基酸对稳定性的贡献（盐桥）
    charged_total = comp['K'] + comp['R'] + comp['D'] + comp['E'] + comp['H']
    if charged_total > 25:
        score += 5

    # 芳香族堆积（疏水核心中）
    if aromatic := comp['F'] + comp['Y'] + comp['W'] > 0:
        score += min(aromatic * 2, 15)

    return max(0, min(100, score))


def analyze_sequence(seq_id, seq):
    """分析单条序列"""
    length = len(seq)
    comp = calculate_composition(seq)

    expressibility = calculate_expressibility_score(seq, comp)
    solubility = calculate_solubility_score(seq, comp)
    structure = calculate_structure_score(seq, comp)
    stability = calculate_stability_score(seq, comp)

    # 综合评分（加权平均）
    # 权重：表达性20%，水溶性25%，结构25%，稳定性30%
    total_score = (expressibility * 0.20 +
                   solubility * 0.25 +
                   structure * 0.25 +
                   stability * 0.30)

    return {
        'id': seq_id,
        'sequence': seq,
        'length': length,
        'expressibility': round(expressibility, 2),
        'solubility': round(solubility, 2),
        'structure': round(structure, 2),
        'stability': round(stability, 2),
        'total_score': round(total_score, 2)
    }


def main():
    # 文件路径
    work_dir = r"c:\Users\lenovo\Desktop\pmpnn"
    fasta_file = os.path.join(work_dir, "chain_a_sequences.fasta")
    output_file = os.path.join(work_dir, "protein_scores.txt")
    top50_file = os.path.join(work_dir, "top50_sequences.fasta")

    # 读取序列
    print("正在读取FASTA文件...")
    sequences = read_fasta(fasta_file)
    print(f"共读取 {len(sequences)} 条序列\n")

    # 分析所有序列
    print("正在分析序列...")
    results = []
    for i, (seq_id, seq) in enumerate(sequences.items(), 1):
        if i % 100 == 0:
            print(f"已分析 {i}/{len(sequences)} 条序列...")
        result = analyze_sequence(seq_id, seq)
        results.append(result)

    # 按总分排序
    results.sort(key=lambda x: x['total_score'], reverse=True)

    # 写入评分表
    print(f"\n正在生成评分表...")
    with open(output_file, 'w', encoding='utf-8') as f:
        # 写入表头
        header = f"{'排名':<6}{'序列ID':<55}{'长度':<8}{'表达性':<10}{'水溶性':<10}{'结构':<10}{'稳定性':<10}{'总分':<10}\n"
        f.write("=" * 120 + "\n")
        f.write("蛋白质综合评估打分表\n")
        f.write("评分说明: 表达性(0-100)、水溶性(0-100)、结构合理性(0-100)、稳定性(0-100)、总分(0-100)\n")
        f.write("权重分配: 表达性20%、水溶性25%、结构25%、稳定性30%\n")
        f.write("=" * 120 + "\n")
        f.write(header)
        f.write("-" * 120 + "\n")

        # 写入所有结果
        for rank, result in enumerate(results, 1):
            line = f"{rank:<6}{result['id']:<55}{result['length']:<8}"
            line += f"{result['expressibility']:<10.2f}{result['solubility']:<10.2f}"
            line += f"{result['structure']:<10.2f}{result['stability']:<10.2f}"
            line += f"{result['total_score']:<10.2f}\n"
            f.write(line)

    print(f"评分表已保存到: {output_file}")

    # 提取前50条序列
    print(f"\n正在提取评分最高的前50条序列...")
    with open(top50_file, 'w', encoding='utf-8') as f:
        for rank, result in enumerate(results[:50], 1):
            f.write(f">{result['id']} [排名:{rank} 总分:{result['total_score']:.2f}]\n")
            seq = result['sequence']
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

    print(f"Top 50序列已保存到: {top50_file}")

    # 打印统计信息
    print("\n" + "=" * 60)
    print("统计信息")
    print("=" * 60)
    scores = [r['total_score'] for r in results]
    print(f"总分范围: {min(scores):.2f} - {max(scores):.2f}")
    print(f"总分平均值: {sum(scores)/len(scores):.2f}")
    print(f"前50名最低分: {results[49]['total_score']:.2f}")

    print("\n前10名序列:")
    print("-" * 100)
    for rank, result in enumerate(results[:10], 1):
        print(f"#{rank:3d} | {result['id']:<45} | "
              f"总分:{result['total_score']:5.2f} | "
              f"表达:{result['expressibility']:5.2f} | "
              f"溶解:{result['solubility']:5.2f} | "
              f"结构:{result['structure']:5.2f} | "
              f"稳定:{result['stability']:5.2f}")


if __name__ == "__main__":
    main()
