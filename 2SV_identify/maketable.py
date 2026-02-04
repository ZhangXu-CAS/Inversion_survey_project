import os
import pandas as pd

base_dir = "/project/6004758/zhangxu/ALLinv_verified"  # 根目录，每个样本一个文件夹
all_data = []

for sample_id in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, sample_id, f"{sample_id}.inv.verify.csv")
    if os.path.exists(sample_path):
        df = pd.read_csv(sample_path)
        df["SampleID"] = sample_id   # 添加样本ID
        all_data.append(df)

# 合并所有样本
master = pd.concat(all_data, ignore_index=True)

# 保存
master.to_csv("Inversion_Master_Table.csv", index=False)
print("Master table saved to Inversion_Master_Table.csv")
