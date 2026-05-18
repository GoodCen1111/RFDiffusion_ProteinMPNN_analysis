"""
蛋白质刚柔性分析脚本
筛选：好表达 + 稳定 + 适度柔性的序列
"""

import os
import re
from collections import OrderedDict

# 氨基酸柔刚性质定义
FLEXIBLE_AA = set('GASP')       # 柔性氨基酸
RIGID_AA = set('FWYCVILM')      # 刚性氨基酸
SEMI_FLEXIBLE = set('NDETQKRH') # 半柔性氨基酸
HYDROPHOBIC_CORE = set('VILMFWY')  # 疏水核心
POLAR_AA = set('NQSTYC')        # 极性氨基酸
CHARGED_AA = set('DEKRH')       # 带电氨基酸
PROLINE = set('P')              # 脯氨酸（破坏二级结构，增加柔性）
GLYCINE = set('G')              # 甘氨酸（高柔性）
CYSTEINE = set('C')             # 半胱氨酸（二硫键）


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
                current_id = line[1:].split('[')[0].strip()
                current_seq = []
            else:
                current_seq.append(line)

        if current_id:
            sequences[current_id] = ''.join(current_seq)

    return sequences


def read_scores(score_file):
    """读取评分文件"""
    scores = {}
    with open(score_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('|') and not line.startswith('|--'):
                parts = [p.strip() for p in line.split('|')]
                if len(parts) >= 9 and parts[2] and not parts[2].startswith('排名'):
                    try:
                        rank = int(parts[1])
                        seq_id = parts[2]
                        scores[seq_id] = {
                            'rank': rank,
                            'expressibility': float(parts[4]),
                            'solubility': float(parts[5]),
                            'structure': float(parts[6]),
                            'stability': float(parts[7]),
                            'total': float(parts[8])
                        }
                    except (ValueError, IndexError):
                        continue
    return scores


def analyze_flexibility_rigidity(seq):
    """
    分析蛋白质刚柔性
    返回详细的柔刚性分析结果
    """
    length = len(seq)
    if length == 0:
        return None

    # 计算各类型氨基酸数量
    flexible_count = sum(1 for aa in seq if aa in FLEXIBLE_AA)
    rigid_count = sum(1 for aa in seq if aa in RIGID_AA)
    semi_flexible_count = sum(1 for aa in seq if aa in SEMI_FLEXIBLE)
    glycine_count = seq.count('G')
    proline_count = seq.count('P')
    cysteine_count = seq.count('C')

    # 比例计算
    flex_ratio = flexible_count / length * 100
    rigid_ratio = rigid_count / length * 100
    semi_ratio = semi_flexible_count / length * 100
    gly_ratio = glycine_count / length * 100
    pro_ratio = proline_count / length * 100
    cys_ratio = cysteine_count / length * 100

    # 疏水核心分析（刚性区域）
    hydrophobic_core = sum(1 for aa in seq if aa in HYDROPHOBIC_CORE)
    core_ratio = hydrophobic_core / length * 100

    # 极性氨基酸
    polar_count = sum(1 for aa in seq if aa in POLAR_AA)
    polar_ratio = polar_count / length * 100

    # 带电氨基酸
    charged_count = sum(1 for aa in seq if aa in CHARGED_AA)
    charged_ratio = charged_count / length * 100

    # 柔性-刚性比
    flex_rigid_ratio = flexible_count / rigid_count if rigid_count > 0 else flexible_count

    # Pro-Gly分析（转角区域，高柔性）
    pg_motifs = len(re.findall(r'P.G|G.P|G.G|P.P', seq))

    # 连续柔性氨基酸
    flex_stretches = re.findall(r'[GASP]{5,}', seq)
    max_flex_stretch = max([len(s) for s in flex_stretches]) if flex_stretches else 0

    # 连续刚性氨基酸
    rigid_stretches = re.findall(r'[FWYCVILM]{5,}', seq)
    max_rigid_stretch = max([len(s) for s in rigid_stretches]) if rigid_stretches else 0

    # 二硫键潜力
    disulfide_potential = cysteine_count // 2

    # 表面暴露vs埋藏倾向
    # 表面倾向：带电、极性、小的疏水
    surface_score = charged_ratio * 2 + polar_ratio + flexible_count / length * 100
    # 核心倾向：大的疏水
    core_score = core_ratio * 3

    return {
        'length': length,
        'flexible_count': flexible_count,
        'rigid_count': rigid_count,
        'flex_ratio': flex_ratio,
        'rigid_ratio': rigid_ratio,
        'semi_ratio': semi_ratio,
        'gly_ratio': gly_ratio,
        'pro_ratio': pro_ratio,
        'cys_ratio': cys_ratio,
        'core_ratio': core_ratio,
        'polar_ratio': polar_ratio,
        'charged_ratio': charged_ratio,
        'flex_rigid_ratio': flex_rigid_ratio,
        'pg_motifs': pg_motifs,
        'max_flex_stretch': max_flex_stretch,
        'max_rigid_stretch': max_rigid_stretch,
        'disulfide_potential': disulfide_potential,
        'surface_score': surface_score,
        'core_score': core_score
    }


def calculate_expression_score(analysis):
    """
    表达性评分（基于刚柔性）
    - 适度的柔性好表达
    - 过多的柔性不利于稳定表达
    - 适当的刚性有助于折叠
    """
    score = 100
    flex_ratio = analysis['flex_ratio']
    rigid_ratio = analysis['rigid_ratio']
    pro_ratio = analysis['pro_ratio']
    gly_ratio = analysis['gly_ratio']

    # 理想柔性比例范围：25-40%
    if flex_ratio < 20:
        score -= (20 - flex_ratio) * 2  # 太刚性，不利于表达
    elif flex_ratio > 45:
        score -= (flex_ratio - 45) * 1.5  # 太柔性，可能导致错误折叠

    # 理想刚性比例范围：20-35%
    if rigid_ratio < 15:
        score -= (15 - rigid_ratio) * 2
    elif rigid_ratio > 40:
        score -= (rigid_ratio - 40) * 1.5

    # 脯氨酸过多不利于表达（翻译停顿）
    if pro_ratio > 8:
        score -= (pro_ratio - 8) * 3

    # 甘氨酸过多降低稳定性
    if gly_ratio > 12:
        score -= (gly_ratio - 12) * 2

    # 适度的柔刚比（0.8-1.5）
    fr_ratio = analysis['flex_rigid_ratio']
    if fr_ratio < 0.6:
        score -= (0.6 - fr_ratio) * 20
    elif fr_ratio > 2.0:
        score -= (fr_ratio - 2.0) * 15

    return max(0, min(100, score))


def calculate_flexibility_score(analysis):
    """
    柔性评分
    - 需要适度柔性用于功能运动
    - 过多柔性导致无序/不稳定
    - 过少柔性导致僵硬/无活性
    """
    score = 100
    flex_ratio = analysis['flex_ratio']
    gly_ratio = analysis['gly_ratio']
    pro_ratio = analysis['pro_ratio']

    # 适度柔性评分（基于功能蛋白的理想范围）
    # 功能蛋白通常有25-40%的柔性氨基酸
    if 25 <= flex_ratio <= 40:
        score += 10
    elif flex_ratio < 20:
        score -= (20 - flex_ratio) * 3
    elif flex_ratio > 50:
        score -= (flex_ratio - 50) * 2

    # 甘氨酸加分（高度柔性）
    if 3 <= gly_ratio <= 10:
        score += 10
    elif gly_ratio > 12:
        score -= (gly_ratio - 12) * 3

    # 脯氨酸加分（转角/铰链区域）
    if 3 <= pro_ratio <= 8:
        score += 10
    elif pro_ratio > 10:
        score -= (pro_ratio - 10) * 2

    # Pro-Gly motif加分（经典柔性转角）
    pg_score = min(analysis['pg_motifs'] * 3, 15)
    score += pg_score

    # 连续柔性序列过多减分
    if analysis['max_flex_stretch'] > 10:
        score -= (analysis['max_flex_stretch'] - 10) * 2

    return max(0, min(100, score))


def calculate_stability_score_fr(analysis):
    """
    稳定性评分（基于刚柔性）
    - 需要足够的刚性核心
    - 二硫键加分
    - 柔刚平衡
    """
    score = 100
    core_ratio = analysis['core_ratio']
    cys_ratio = analysis['cys_ratio']
    rigid_ratio = analysis['rigid_ratio']
    disulfide = analysis['disulfide_potential']

    # 疏水核心评分（理想范围20-35%）
    if 20 <= core_ratio <= 35:
        score += 15
    elif core_ratio < 15:
        score -= (15 - core_ratio) * 3
    elif core_ratio > 40:
        score -= (core_ratio - 40) * 2

    # 刚性比例加分（稳定性）
    if 20 <= rigid_ratio <= 35:
        score += 10

    # 二硫键评分
    if disulfide >= 2:
        score += min(disulfide * 5, 20)
    elif disulfide == 1:
        score += 5

    # 连续的刚性序列加分
    if analysis['max_rigid_stretch'] >= 5:
        score += min(analysis['max_rigid_stretch'], 15)

    # 过多的连续柔性减分
    if analysis['max_flex_stretch'] > 12:
        score -= (analysis['max_flex_stretch'] - 12) * 3

    # 带电残基（盐桥）加分
    charged = analysis['charged_ratio']
    if charged >= 20:
        score += min((charged - 20) * 0.5, 10)

    return max(0, min(100, score))


def calculate_balance_score(analysis):
    """
    柔刚平衡评分
    好表达 + 稳定 + 适度柔性 = 最佳平衡
    """
    score = 100

    flex_ratio = analysis['flex_ratio']
    rigid_ratio = analysis['rigid_ratio']
    flex_rigid_ratio = analysis['flex_rigid_ratio']

    # 理想柔刚比例（约1:1）
    if 0.7 <= flex_rigid_ratio <= 1.5:
        score += 15
    elif 0.5 <= flex_rigid_ratio < 0.7:
        score += 10
    elif 1.5 < flex_rigid_ratio <= 2.0:
        score += 10
    else:
        score -= 10

    # 表面/核心平衡
    surface = analysis['surface_score']
    core = analysis['core_score']
    if surface > 0 and core > 0:
        ratio = surface / core
        if 0.5 <= ratio <= 2.0:
            score += 10

    # 二硫键稳定性贡献
    disulfide = analysis['disulfide_potential']
    if disulfide >= 2:
        score += 15
    elif disulfide == 1:
        score += 5

    return max(0, min(100, score))


def analyze_sequence(seq_id, seq, existing_scores):
    """综合分析单条序列"""
    analysis = analyze_flexibility_rigidity(seq)

    expression_score = calculate_expression_score(analysis)
    flexibility_score = calculate_flexibility_score(analysis)
    stability_score = calculate_stability_score_fr(analysis)
    balance_score = calculate_balance_score(analysis)

    # 合并之前的评分
    prev_scores = existing_scores.get(seq_id, {})

    # 综合评分（考虑之前的表达性、水溶性、结构、稳定性 + 新的柔刚分析）
    total_score = (
        expression_score * 0.25 +
        flexibility_score * 0.20 +
        stability_score * 0.25 +
        balance_score * 0.15 +
        prev_scores.get('expressibility', 50) * 0.05 +
        prev_scores.get('solubility', 50) * 0.05 +
        prev_scores.get('stability', 50) * 0.05
    )

    return {
        'id': seq_id,
        'sequence': seq,
        'length': analysis['length'],
        'analysis': analysis,
        'expression_score': round(expression_score, 2),
        'flexibility_score': round(flexibility_score, 2),
        'stability_score': round(stability_score, 2),
        'balance_score': round(balance_score, 2),
        'total_score': round(total_score, 2),
        'prev_total': prev_scores.get('total', 0)
    }


def main():
    work_dir = r"c:\Users\lenovo\Desktop\pmpnn"
    fasta_file = os.path.join(work_dir, "chain_a_sequences.fasta")
    score_file = os.path.join(work_dir, "protein_scores.txt")
    output_file = os.path.join(work_dir, "flexibility_analysis.txt")
    top50_flex_file = os.path.join(work_dir, "top50_flexible_sequences.fasta")

    # 读取数据
    print("正在读取数据...")
    sequences = read_fasta(fasta_file)
    existing_scores = read_scores(score_file)
    print(f"共读取 {len(sequences)} 条序列\n")

    # 分析所有序列
    print("正在分析刚柔性...")
    results = []
    for i, (seq_id, seq) in enumerate(sequences.items(), 1):
        if i % 100 == 0:
            print(f"已分析 {i}/{len(sequences)} 条序列...")
        result = analyze_sequence(seq_id, seq, existing_scores)
        results.append(result)

    # 按总分排序
    results.sort(key=lambda x: x['total_score'], reverse=True)

    # 生成详细分析报告
    print(f"\n正在生成分析报告...")
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("=" * 140 + "\n")
        f.write("蛋白质刚柔性综合分析报告\n")
        f.write("=" * 140 + "\n\n")

        f.write("【评分说明】\n")
        f.write("1. 表达性评分(0-100): 基于柔刚比例、Pro/Gly含量、柔刚平衡\n")
        f.write("2. 柔性评分(0-100): 适度柔性好功能，过多/过少都不好\n")
        f.write("3. 稳定性评分(0-100): 疏水核心、二硫键、柔刚平衡\n")
        f.write("4. 平衡评分(0-100): 表面/核心平衡、二硫键稳定性\n")
        f.write("5. 综合总分: 综合上述所有指标\n\n")

        f.write("【理想参数范围】\n")
        f.write("- 柔性氨基酸比例: 25-40%\n")
        f.write("- 刚性氨基酸比例: 20-35%\n")
        f.write("- 疏水核心比例: 20-35%\n")
        f.write("- 柔刚比: 0.7-1.5\n")
        f.write("- 甘氨酸比例: 3-10%\n")
        f.write("- 脯氨酸比例: 3-8%\n")
        f.write("- 二硫键数量: >=2\n\n")

        f.write("=" * 140 + "\n")
        header = f"{'排名':<5}{'序列ID':<50}{'长度':<6}"
        header += f"{'表达性':<8}{'柔性':<8}{'稳定性':<8}{'平衡':<8}{'总分':<8}{'之前总分':<10}\n"
        f.write(header)
        f.write("-" * 140 + "\n")

        for rank, r in enumerate(results, 1):
            a = r['analysis']
            line = f"{rank:<5}{r['id']:<50}{r['length']:<6}"
            line += f"{r['expression_score']:<8.2f}{r['flexibility_score']:<8.2f}"
            line += f"{r['stability_score']:<8.2f}{r['balance_score']:<8.2f}"
            line += f"{r['total_score']:<8.2f}{r['prev_total']:<10.2f}\n"
            f.write(line)

        # 添加详细分析部分
        f.write("\n" + "=" * 140 + "\n")
        f.write("前50名序列详细刚柔性分析\n")
        f.write("=" * 140 + "\n\n")

        for rank, r in enumerate(results[:50], 1):
            a = r['analysis']
            f.write(f"#{rank} {r['id']}\n")
            f.write(f"   总分: {r['total_score']:.2f} | 长度: {r['length']} aa\n")
            f.write(f"   表达性: {r['expression_score']:.2f} | 柔性: {r['flexibility_score']:.2f} | 稳定性: {r['stability_score']:.2f} | 平衡: {r['balance_score']:.2f}\n")
            f.write(f"   柔刚比: {a['flex_ratio']:.1f}% / {a['rigid_ratio']:.1f}% = {a['flex_rigid_ratio']:.2f}\n")
            f.write(f"   Gly: {a['gly_ratio']:.1f}% | Pro: {a['pro_ratio']:.1f}% | Cys: {a['cys_ratio']:.1f}%\n")
            f.write(f"   疏水核心: {a['core_ratio']:.1f}% | 极性: {a['polar_ratio']:.1f}% | 带电: {a['charged_ratio']:.1f}%\n")
            f.write(f"   二硫键潜力: {a['disulfide_potential']} | Pro-Gly motifs: {a['pg_motifs']}\n")
            f.write(f"   最大连续柔性: {a['max_flex_stretch']} | 最大连续刚性: {a['max_rigid_stretch']}\n")
            f.write(f"   序列: {r['sequence'][:60]}...\n\n")

    print(f"分析报告已保存到: {output_file}")

    # 提取前50条柔性好的序列
    print(f"\n正在提取前50条最优序列...")
    with open(top50_flex_file, 'w', encoding='utf-8') as f:
        for rank, r in enumerate(results[:50], 1):
            a = r['analysis']
            info = (f"排名:{rank} 总分:{r['total_score']:.2f} "
                   f"表达:{r['expression_score']:.1f} "
                   f"柔性:{r['flexibility_score']:.1f} "
                   f"稳定:{r['stability_score']:.1f} "
                   f"柔刚比:{a['flex_rigid_ratio']:.2f} "
                   f"二硫键:{a['disulfide_potential']}")
            f.write(f">{r['id']} [{info}]\n")
            seq = r['sequence']
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

    print(f"Top 50序列已保存到: {top50_flex_file}")

    # 统计信息
    print("\n" + "=" * 80)
    print("分析结果统计")
    print("=" * 80)
    scores = [r['total_score'] for r in results]
    print(f"综合总分范围: {min(scores):.2f} - {max(scores):.2f}")
    print(f"综合总分平均值: {sum(scores)/len(scores):.2f}")
    print(f"Top 50门槛分: {results[49]['total_score']:.2f}")

    # Top 10详细信息
    print("\n" + "=" * 100)
    print("Top 10 柔性好、稳定、易表达的序列")
    print("=" * 100)
    print(f"{'排名':<4} {'序列ID':<45} {'总分':<8} {'表达性':<8} {'柔性':<8} {'稳定性':<8} {'平衡':<8} {'柔刚比':<8} {'二硫键':<6}")
    print("-" * 100)

    for rank, r in enumerate(results[:10], 1):
        a = r['analysis']
        print(f"{rank:<4} {r['id']:<45} {r['total_score']:<8.2f} "
              f"{r['expression_score']:<8.2f} {r['flexibility_score']:<8.2f} "
              f"{r['stability_score']:<8.2f} {r['balance_score']:<8.2f} "
              f"{a['flex_rigid_ratio']:<8.2f} {a['disulfide_potential']:<6}")


if __name__ == "__main__":
    main()
