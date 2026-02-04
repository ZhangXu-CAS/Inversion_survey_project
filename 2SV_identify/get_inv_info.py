import os
import sys
import pandas as pd

def parse_inv(inv_file):
    inversions = []
    with open(inv_file) as f:
        for line in f:
            if line.startswith("#"):
                parts = line.strip().split()
                scaffold = parts[1]
                ref_start = int(parts[2])
                ref_end = int(parts[3])
                qry_start = int(parts[6])
                qry_end = int(parts[7])
                ref_length = abs(ref_end - ref_start) + 1
                inversions.append([scaffold, ref_start, ref_end, ref_length, qry_start, qry_end])
    return inversions

if __name__ == "__main__":
    sample = sys.argv[1]
    inv_file = os.path.join(sample, "syri_results_reminimap", "invOut.txt")
    out_txt = os.path.join(sample, f"{sample}_inversions.bed")

    # 解析 inversion 信息
    inv_list = parse_inv(inv_file)
    # 在每行前面加上 sample ID
    inv_list = [[sample] + inv for inv in inv_list]

    # 创建 DataFrame
    df = pd.DataFrame(inv_list, columns=["sample", "scaffold", "hap1_start", "hap1_end", "inv_length", "hap2_start", "hap2_end"])

    # 保存表格
    df.to_csv(out_txt, sep="\t", index=False)