import os
os.chdir("E:/experiment/data analysis/pY of other species/Ecoli pY")
# 在 Jupyter Notebook 中直接运行这部分代码

def read_fasta(file_path):
    """读取FASTA文件"""
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_id is not None:
            sequences[current_id] = ''.join(current_seq)
    
    print(f"成功读取 {len(sequences)} 条蛋白序列")
    return sequences

def extract_peptides(sequences, center_aa, window_size=6):
    """提取以特定氨基酸为中心的肽段"""
    peptides = set()
    total_sites = 0
    extracted_sites = 0
    
    center_aa = center_aa.upper()
    
    for prot_id, seq in sequences.items():
        seq_upper = seq.upper()
        
        for i in range(len(seq_upper)):
            if seq_upper[i] == center_aa:
                total_sites += 1
                
                start_pos = i - window_size
                end_pos = i + window_size + 1
                
                if start_pos >= 0 and end_pos <= len(seq_upper):
                    peptide = seq_upper[start_pos:end_pos]
                    if len(peptide) == 2 * window_size + 1:
                        peptides.add(peptide)
                        extracted_sites += 1
    
    print(f"中心氨基酸 {center_aa}:")
    print(f"  总位点数: {total_sites}")
    print(f"  成功提取位点: {extracted_sites}")
    print(f"  去重后肽段数: {len(peptides)}")
    
    return peptides

def save_peptides(peptides, output_file, center_aa, window_size):
    """保存肽段到文件"""
    with open(output_file, 'w') as f:
        f.write(f"# 中心氨基酸: {center_aa}\n")
        f.write(f"# 窗口大小: ±{window_size}\n")
        f.write(f"# 肽段长度: {2*window_size + 1}\n")
        f.write(f"# 肽段数量: {len(peptides)}\n")
        f.write("#" * 60 + "\n")
        
        sorted_peptides = sorted(peptides)
        for i, peptide in enumerate(sorted_peptides, 1):
            marked_peptide = list(peptide)
            marked_peptide[window_size] = marked_peptide[window_size].lower()
            marked_peptide = ''.join(marked_peptide)
            
            f.write(f"{i:6d}\t{peptide}\t{marked_peptide}\n")
    
    print(f"肽段已保存到: {output_file}")

# === 在这里设置你的参数 ===
fasta_file = "ecoli_proteins.fasta"  # 替换为你的FASTA文件名
window_size = 6
output_prefix = "peptides"
# ==========================

print("开始处理...")
print(f"输入文件: {fasta_file}")
print(f"窗口大小: ±{window_size}")

# 读取FASTA文件
sequences = read_fasta(fasta_file)

if not sequences:
    print("错误: 没有读取到任何蛋白序列！")
else:
    # 分别提取Y为中心的肽段
    center_amino_acids = ['Y']
    
    for aa in center_amino_acids:
        print(f"\n提取以 {aa} 为中心的肽段...")
        peptides = extract_peptides(sequences, aa, window_size)
        
        # 保存结果
        output_file = f"{output_prefix}_{aa}.txt"
        save_peptides(peptides, output_file, aa, window_size)
    
    print("\n任务完成！")