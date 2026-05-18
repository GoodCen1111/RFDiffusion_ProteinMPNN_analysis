import os
from pathlib import Path

# ===================== 配置项（可修改） =====================
INPUT_FOLDER = "."  # 当前文件夹（无需改动）
OUTPUT_FILE = "A链序列合集.txt"  # 输出文件名
# ============================================================

def extract_chain_a_sequence(pdb_file):
    """
    从PDB文件中提取A链的氨基酸序列
    :param pdb_file: PDB文件路径
    :return: (文件名, 序列) 元组，提取失败返回 None
    """
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
        'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
        'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
        'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }

    sequence = []
    prev_res_num = None  # 去重重复残基

    try:
        with open(pdb_file, 'r', encoding='utf-8') as f:
            for line in f:
                # 只读取ATOM行（标准PDB格式）
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    chain_id = line[21].strip()  # 链ID位置
                    if chain_id != "A":
                        continue

                    res_name = line[17:20].strip()  # 氨基酸三字母码
                    res_num = line[22:26].strip()   # 残基编号

                    # 同一个残基只提取一次
                    if res_num != prev_res_num:
                        if res_name in three_to_one:
                            sequence.append(three_to_one[res_name])
                        prev_res_num = res_num

        seq_str = ''.join(sequence)
        if not seq_str:
            return None
        return Path(pdb_file).stem, seq_str  # 返回文件名（无后缀）+序列

    except Exception as e:
        print(f"❌ 读取失败 {pdb_file.name}: {str(e)}")
        return None

def main():
    print("🔍 开始扫描当前文件夹的PDB文件...")
    pdb_files = list(Path(INPUT_FOLDER).glob("*.pdb"))

    if not pdb_files:
        print("❌ 未找到任何PDB文件！")
        return

    # 存储所有序列
    all_sequences = []
    success_count = 0

    for pdb in pdb_files:
        result = extract_chain_a_sequence(pdb)
        if result:
            name, seq = result
            # FASTA格式
            fasta_entry = f">{name}\n{seq}\n"
            all_sequences.append(fasta_entry)
            success_count += 1
            print(f"✅ 提取成功: {pdb.name}")
        else:
            print(f"⚠️  无A链或提取失败: {pdb.name}")

    # 写入文件
    with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
        f.writelines(all_sequences)

    print(f"\n🎉 任务完成！")
    print(f"📊 总计处理: {len(pdb_files)} 个PDB文件")
    print(f"✅ 成功提取: {success_count} 条A链序列")
    print(f"📄 序列已保存到: {OUTPUT_FILE}")

if __name__ == "__main__":
    main()