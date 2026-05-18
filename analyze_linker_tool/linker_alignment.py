"""
提取Linker序列并进行多序列比对可视化
"""

import os
from collections import OrderedDict

# 可视化库
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    HAS_PLOT = True
except ImportError:
    HAS_PLOT = False

# ============ 序列定义 ============
# 对照序列 (cIgg2 linker) - 用户提供
CONTROL_SEQUENCE = {
    'name': 'cIgg2 (Control)',
    'linker': 'EPKIPQPQPKPQPQPQPQPKPQPKPEPECTCPKCP'
}

# 用户序列
USER_SEQUENCE = {
    'name': 'User Sequence (#104)',
    'linker': 'IVKSIKKLLSKIKNLKELKKKKKKKKKKAKGGKTLPVSNPGDPVG'
}

# 固定端序列用于验证
BJ_SEQ = "KKNGYPLDRNGKTTECSGVNAIAPHYCNSECTKVYYAESGYCCWGACYCFGLEDDKPIGPMKDITKKYCDVQ"
LqhIT_SEQ = "DAYIAKNYNCVYECFRDAYCNELCTKNGASSGYCQWAGKYGNACWCYALPDNVPIRVPGKCHR"


def read_fasta(fasta_file):
    """读取FASTA文件"""
    sequences = OrderedDict()
    current_id = None
    current_seq = []

    with open(fasta_file, 'r', encoding='utf-8') as f:
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


def extract_linker(full_seq):
    """提取linker序列

    序列结构: [BJ] + [Linker] + [LqhαIT]
    BJ结束: YCDVQ (完整序列的56个氨基酸)
    IT开始: DAYIAKNYNC...
    """
    # BJ序列结束位置 (56个氨基酸)
    BJ_LENGTH = 56
    # 或者用模式匹配
    bj_pattern = "KKNGYPLDRNGKTTECSGVNAIAPHYCNSECTKVYYAESGYCCWGACYCFGLEDDKPIGPMKDITKKYCDVQ"

    # IT序列开始位置 (结尾是CHR不是CHRK)
    it_pattern = "DAYIAKNYNCVYECFRDAYCNELCTKNGASSGYCQWAGKYGNACWCYALPDNVPIRVPGKCHR"

    bj_idx = full_seq.find(bj_pattern)
    it_idx = full_seq.find(it_pattern)

    if bj_idx == -1:
        print(f"  WARNING: BJ pattern not found!")
        return None
    if it_idx == -1:
        print(f"  WARNING: IT pattern not found!")
        return None

    # linker从BJ结束后开始，到IT前结束
    linker = full_seq[bj_idx + len(bj_pattern):it_idx]

    return linker


def validate_extraction():
    """验证序列提取是否正确"""
    print("\n" + "=" * 60)
    print("Validating Linker Extraction")
    print("=" * 60)

    # 读取序列
    fasta_file = r"c:\Users\lenovo\Desktop\pmpnn\top50_sequences.fasta"
    sequences = read_fasta(fasta_file)

    print(f"\nTotal sequences loaded: {len(sequences)}")

    # 检查前3条
    for i, (seq_id, seq) in enumerate(list(sequences.items())[:3]):
        print(f"\n--- Sequence #{i+1}: {seq_id[:50]} ---")
        print(f"Total length: {len(seq)} aa")

        # 验证BJ序列
        bj_found = seq.find(BJ_SEQ[:30])
        print(f"BJ sequence found at: {bj_found}")

        # 验证IT序列
        it_found = seq.find(LqhIT_SEQ[:30])
        print(f"LqhαIT sequence found at: {it_found}")

        # 提取linker
        linker = extract_linker(seq)
        if linker:
            print(f"Linker extracted: {len(linker)} aa")
            print(f"Linker sequence: {linker}")
        else:
            print("WARNING: Failed to extract linker!")


def get_consensus_sequence(sequences):
    """获取共有序列"""
    if not sequences:
        return ""

    min_len = min(len(seq) for seq in sequences)
    consensus = []

    for i in range(min_len):
        counts = {}
        for seq in sequences:
            aa = seq[i]
            counts[aa] = counts.get(aa, 0) + 1
        consensus_aa = max(counts.items(), key=lambda x: x[1])[0]
        consensus.append(consensus_aa)

    return ''.join(consensus)


def plot_alignment(all_linkers, seq_ids, output_file):
    """绘制多序列比对图"""
    if not HAS_PLOT:
        print("需要 matplotlib 来生成比对图")
        return

    print("\nGenerating linker sequence alignment...")

    # 设置Arial字体
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['axes.unicode_minus'] = False

    n_seqs = len(all_linkers)
    if n_seqs == 0:
        print("No linkers extracted!")
        return

    max_len = max(len(linker) for linker in all_linkers)

    # 颜色映射
    aa_colors = {
        'K': '#FF6B6B', 'R': '#FF6B6B',  # 碱性
        'D': '#4ECDC4', 'E': '#4ECDC4',  # 酸性
        'A': '#95E1D3', 'V': '#95E1D3', 'L': '#95E1D3', 'I': '#95E1D3', 'M': '#95E1D3',  # 疏水
        'F': '#A8E6CF', 'W': '#A8E6CF', 'Y': '#A8E6CF',  # 芳香族
        'S': '#FFE66D', 'T': '#FFE66D',  # 极性
        'N': '#DDA0DD', 'Q': '#DDA0DD',  #酰胺
        'G': '#E0E0E0',  # 甘氨酸
        'P': '#C0C0C0',  # 脯氨酸
        'C': '#FFD700',  # 半胱氨酸
    }

    # 图例元素
    legend_elements = [
        patches.Patch(facecolor='#FF6B6B', label='Basic (K,R)'),
        patches.Patch(facecolor='#4ECDC4', label='Acidic (D,E)'),
        patches.Patch(facecolor='#95E1D3', label='Hydrophobic (A,V,L,I,M)'),
        patches.Patch(facecolor='#A8E6CF', label='Aromatic (F,W,Y)'),
        patches.Patch(facecolor='#FFE66D', label='Polar (S,T)'),
        patches.Patch(facecolor='#DDA0DD', label='Amide (N,Q)'),
        patches.Patch(facecolor='#E0E0E0', label='Glycine (G)'),
        patches.Patch(facecolor='#C0C0C0', label='Proline (P)'),
        patches.Patch(facecolor='#FFD700', label='Cysteine (C)'),
    ]

    # 创建图形
    fig = plt.figure(figsize=(24, 14))

    # 左侧区域 - 图例
    ax_legend = fig.add_axes([0.02, 0.12, 0.18, 0.78])
    ax_legend.axis('off')
    ax_legend.legend(handles=legend_elements, loc='center', ncol=1,
                    prop={'family': 'Arial', 'size': 10}, frameon=True,
                    title='Amino Acid Color Code', title_fontsize=11)
    ax_legend.set_xlim(0, 1)
    ax_legend.set_ylim(0, 1)

    # 中间区域 - 序列比对
    ax_seq = fig.add_axes([0.22, 0.12, 0.76, 0.78])

    cell_width = 1.0
    cell_height = 1.0
    row_spacing = 0.15

    total_height = n_seqs * (cell_height + row_spacing) + 2
    total_width = max_len + 30

    # 绘制每个序列
    for i, (linker, seq_id) in enumerate(zip(all_linkers, seq_ids)):
        y_pos = total_height - i * (cell_height + row_spacing) - cell_height

        # 序列ID显示
        ax_seq.text(-1, y_pos + cell_height/2, f"#{i+1}", ha='right', va='center',
                   fontsize=9, fontfamily='Arial', fontweight='bold')

        # 根据是否为对照序列设置标签
        if 'cIgg2' in seq_id:
            label_text = seq_id
            color = '#E74C3C'  # 红色标注对照
        elif 'User' in seq_id:
            label_text = seq_id
            color = '#3498DB'  # 蓝色标注用户序列
        else:
            label_text = seq_id[:22] + '...' if len(seq_id) > 22 else seq_id
            color = 'black'

        ax_seq.text(0, y_pos + cell_height/2, label_text, ha='left', va='center',
                   fontsize=7, fontfamily='Arial', style='italic', color=color)

        # 绘制序列格子
        for j, aa in enumerate(linker):
            color = aa_colors.get(aa, '#FFFFFF')
            rect = patches.Rectangle((j + 30, y_pos), cell_width, cell_height,
                                     facecolor=color, edgecolor='#888888', linewidth=0.3)
            ax_seq.add_patch(rect)
            ax_seq.text(j + 30.5, y_pos + cell_height/2, aa, ha='center', va='center',
                       fontsize=8, fontweight='bold', fontfamily='Arial')

        # 灰色填充空白
        for j in range(len(linker), max_len):
            rect = patches.Rectangle((j + 30, y_pos), cell_width, cell_height,
                                     facecolor='#F8F8F8', edgecolor='#DDDDDD', linewidth=0.2)
            ax_seq.add_patch(rect)

    # 绘制共有序列
    consensus = get_consensus_sequence([l for l, s in zip(all_linkers, seq_ids) if 'cIgg2' not in s and 'User' not in s])
    y_pos_consensus = total_height - n_seqs * (cell_height + row_spacing) - 2

    ax_seq.text(-1, y_pos_consensus + cell_height/2, "Con.", ha='right', va='center',
               fontsize=9, fontweight='bold', fontfamily='Arial', color='#333333')
    ax_seq.text(0, y_pos_consensus + cell_height/2, "Consensus (Top 10)", ha='left', va='center',
               fontsize=8, fontweight='bold', fontfamily='Arial', color='#333333')

    for j, aa in enumerate(consensus):
        color = aa_colors.get(aa, '#FFFFFF')
        rect = patches.Rectangle((j + 30, y_pos_consensus), cell_width, cell_height,
                                 facecolor=color, edgecolor='#333333', linewidth=0.8)
        ax_seq.add_patch(rect)
        ax_seq.text(j + 30.5, y_pos_consensus + cell_height/2, aa, ha='center', va='center',
                   fontsize=9, fontweight='bold', fontfamily='Arial', color='black')

    # 设置坐标轴
    ax_seq.set_xlim(-2, max_len + 32)
    ax_seq.set_ylim(y_pos_consensus - 1, total_height + 1)
    ax_seq.set_aspect('equal')
    ax_seq.axis('off')

    # 刻度线
    for pos in range(0, max_len, 10):
        ax_seq.axvline(x=pos + 30, color='lightgray', linestyle='-', linewidth=0.3, zorder=0)

    # 标题
    fig.suptitle(f'Linker Sequence Alignment (Top 10 + User + Control: cIgg2)',
                fontsize=16, fontweight='bold', fontfamily='Arial', y=0.95)

    # 底部统计
    top10_linkers = [l for l, s in zip(all_linkers, seq_ids) if 'cIgg2' not in s and 'User' not in s]
    lengths = [len(linker) for linker in top10_linkers]
    user_len = len([l for l, s in zip(all_linkers, seq_ids) if 'User' in s][0])
    cIgg2_len = len([l for l, s in zip(all_linkers, seq_ids) if 'cIgg2' in s][0])

    stats_text = (f"Top 10: Length range = {min(lengths)}-{max(lengths)} aa (avg={sum(lengths)/len(lengths):.1f}) | "
                  f"User (#104): {user_len} aa | "
                  f"cIgg2: {cIgg2_len} aa | "
                  f"Consensus length = {len(consensus)} aa")
    fig.text(0.5, 0.06, stats_text, ha='center', fontsize=10, fontfamily='Arial', style='italic')

    # 保存
    plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"  - Saved: {output_file}")

    return consensus


def main():
    work_dir = r"c:\Users\lenovo\Desktop\pmpnn"
    fasta_file = os.path.join(work_dir, "top50_sequences.fasta")
    output_file = os.path.join(work_dir, "linker_alignment.png")

    print("=" * 60)
    print("Linker Sequence Extraction and Alignment")
    print("=" * 60)

    # 验证提取
    validate_extraction()

    # 读取序列
    print(f"\nReading sequences from: {fasta_file}")
    sequences = read_fasta(fasta_file)
    print(f"Loaded {len(sequences)} sequences")

    # 提取前10条的linker
    print("\nExtracting linkers from top 10 sequences...")
    linkers = []
    seq_ids = []

    for i, (seq_id, seq) in enumerate(list(sequences.items())[:10]):
        linker = extract_linker(seq)
        if linker:
            linkers.append(linker)
            seq_ids.append(seq_id)
            print(f"  #{i+1}: Linker ({len(linker)} aa) - {linker}")
        else:
            print(f"  #{i+1}: FAILED to extract linker!")

    # 添加用户序列
    print(f"\nAdding user sequence (Rank #104)...")
    linkers.append(USER_SEQUENCE['linker'])
    seq_ids.append(USER_SEQUENCE['name'])
    print(f"  User: Linker ({len(USER_SEQUENCE['linker'])} aa) - {USER_SEQUENCE['linker']}")

    # 添加cIgg2对照
    print(f"\nAdding control sequence: cIgg2")
    linkers.append(CONTROL_SEQUENCE['linker'])
    seq_ids.append(CONTROL_SEQUENCE['name'])
    print(f"  cIgg2: Linker ({len(CONTROL_SEQUENCE['linker'])} aa) - {CONTROL_SEQUENCE['linker']}")

    # 生成比对图
    if linkers:
        consensus = plot_alignment(linkers, seq_ids, output_file)

        # 保存序列文件
        linker_file = os.path.join(work_dir, "linker_alignment_sequences.fasta")
        with open(linker_file, 'w', encoding='utf-8') as f:
            f.write(f"# Linker Sequences (Top 10 + User + cIgg2 Control)\n")
            f.write(f"# Consensus: {consensus}\n\n")
            for i, (linker, seq_id) in enumerate(zip(linkers, seq_ids)):
                f.write(f">{seq_id}\n")
                for j in range(0, len(linker), 60):
                    f.write(linker[j:j+60] + "\n")
        print(f"  - Saved: {linker_file}")

        # 统计
        print("\n" + "=" * 60)
        print("Linker Sequence Statistics")
        print("=" * 60)
        top10_linkers = linkers[:-2]  # 排除cIgg2和用户序列
        lengths = [len(l) for l in top10_linkers]
        print(f"Top 10 length range: {min(lengths)} - {max(lengths)} aa")
        print(f"Top 10 average: {sum(lengths)/len(lengths):.1f} aa")
        print(f"User length: {len(linkers[-2])} aa")
        print(f"cIgg2 length: {len(linkers[-1])} aa")
        print(f"\nConsensus sequence: {consensus}")

    print("\n" + "=" * 60)
    print("Done!")
    print("=" * 60)


if __name__ == "__main__":
    main()
