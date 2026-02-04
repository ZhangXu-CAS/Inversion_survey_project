#!/usr/bin/env python3

import sys
import os
import subprocess
from Bio import SeqIO
import pandas as pd
from concurrent.futures import ThreadPoolExecutor

# ========================
# 获取命令行参数
# ========================
if len(sys.argv) != 2:
    print("Usage: python extract_chr.py <sample_id>")
    sys.exit(1)

SAMPLE_ID = sys.argv[1]

# 输入文件
HAP1_FILE = f"{SAMPLE_ID}.hap1_scaffolds_final.fa"
HAP2_FILE = f"{SAMPLE_ID}.hap2_scaffolds_final.fa"

# 输出文件
CHR_NUM_FILE = f"{SAMPLE_ID}_chrn.txt"
HAP1_OUT = f"{SAMPLE_ID}_hap1_chr.fa"
HAP2_OUT = f"{SAMPLE_ID}_hap2_chr.fa"
LOG_FILE = f"{SAMPLE_ID}_extract_chr.log"

# ========================
# Step 1: 提取 scaffold 长度
# ========================
print(f"Processing sample: {SAMPLE_ID}")
# 生成索引
subprocess.run(f"samtools faidx {HAP1_FILE}", shell=True)
subprocess.run(f"samtools faidx {HAP2_FILE}", shell=True)

# 使用生成器避免一次性加载所有数据
def get_scaffold_lengths(file):
    scaffold_lengths = []
    for record in SeqIO.parse(file, "fasta"):
        scaffold_lengths.append((record.id, len(record.seq)))
    return scaffold_lengths

# 获取 hap1 和 hap2 的 scaffold 长度信息
with ThreadPoolExecutor() as executor:
    hap1_lengths_future = executor.submit(get_scaffold_lengths, HAP1_FILE)
    hap2_lengths_future = executor.submit(get_scaffold_lengths, HAP2_FILE)
    
    hap1_lengths = hap1_lengths_future.result()
    hap2_lengths = hap2_lengths_future.result()

# 将数据转换为 DataFrame 并按长度降序排序
hap1_df = pd.DataFrame(hap1_lengths, columns=["scaffold", "length"]).sort_values(by="length", ascending=False).reset_index(drop=True)
hap2_df = pd.DataFrame(hap2_lengths, columns=["scaffold", "length"]).sort_values(by="length", ascending=False).reset_index(drop=True)

# 计算总长度
hap1_total_len = hap1_df["length"].sum()
hap2_total_len = hap2_df["length"].sum()

# ========================
# Step 2: 确定染色体数目
# ========================
# 方法 1: 长度骤减法
chr_num = 1
prev_len = hap1_df.loc[0, "length"]
threshold = 2.0  # 设定阈值，当相邻两个 scaffold 的长度比值大于这个阈值时，认为它们属于不同的染色体

for i in range(1, len(hap1_df)):
    cur_len = hap1_df.loc[i, "length"]
    ratio = prev_len / cur_len
    if ratio > threshold:  # 长度骤减判断
        chr_num = i
        break
    prev_len = cur_len

# 方法 2: 覆盖 95% 长度
if chr_num == 1 or chr_num > 50:  # 如果未检测到骤减，或 scaffold 过多，使用 95% 覆盖法
    sum_len = 0
    for i in range(len(hap1_df)):
        sum_len += hap1_df.loc[i, "length"]
        if sum_len >= 0.95 * hap1_total_len:
            chr_num = i + 1
            break

# 限制最大 scaffold 数量为 50
if chr_num > 50:
    chr_num = 50

# 写入染色体数目
with open(CHR_NUM_FILE, "w") as f:
    f.write(f"{chr_num}\n")

print(f"Detected chromosome number: {chr_num}")

# ========================
# Step 3: 提取染色体 scaffold
# ========================
print("Extracting chromosome-level sequences...")

# 选择前 chr_num 个 scaffold
top_hap1_scaffolds = hap1_df.iloc[:chr_num]["scaffold"].tolist()
top_hap2_scaffolds = hap2_df.iloc[:chr_num]["scaffold"].tolist()

# 提取 hap1 和 hap2 的染色体 scaffold
def extract_scaffolds(file, scaffold_list):
    return [record for record in SeqIO.parse(file, "fasta") if record.id in scaffold_list]

with ThreadPoolExecutor() as executor:
    hap1_chr_records_future = executor.submit(extract_scaffolds, HAP1_FILE, top_hap1_scaffolds)
    hap2_chr_records_future = executor.submit(extract_scaffolds, HAP2_FILE, top_hap2_scaffolds)
    
    hap1_chr_records = hap1_chr_records_future.result()
    hap2_chr_records = hap2_chr_records_future.result()

# 写入染色体序列
SeqIO.write(hap1_chr_records, HAP1_OUT, "fasta")
SeqIO.write(hap2_chr_records, HAP2_OUT, "fasta")

# 生成索引
subprocess.run(f"samtools faidx {HAP1_OUT}", shell=True)
subprocess.run(f"samtools faidx {HAP2_OUT}", shell=True)

# ========================
# Step 4: 计算覆盖度并写入日志
# ========================
print("Calculating coverage and generating log...")

def get_total_length(records):
    """计算提取序列的总长度"""
    return sum(len(record.seq) for record in records)

# 计算提取的染色体序列长度
hap1_chr_len = get_total_length(hap1_chr_records)
hap2_chr_len = get_total_length(hap2_chr_records)

# 计算覆盖率
hap1_cov_rate = round(hap1_chr_len * 100 / hap1_total_len, 2)
hap2_cov_rate = round(hap2_chr_len * 100 / hap2_total_len, 2)

# 生成日志
with open(LOG_FILE, "w") as f:
    f.write(f"{SAMPLE_ID} Chromosome number: {chr_num}\n")
    f.write(f"Hap1 coverage: {hap1_chr_len} / {hap1_total_len} ({hap1_cov_rate}%)\n")
    f.write(f"Hap2 coverage: {hap2_chr_len} / {hap2_total_len} ({hap2_cov_rate}%)\n")

print(f"Log written to {LOG_FILE}")
print("All steps completed successfully!")
