# Inversion Survey Project（中文说明）

[English README](README.md)

本说明由 `scripts/generate_readme.py` 自动生成。

本仓库提供基于单倍型分辨组装的倒位识别流程，适用于多类真核生物数据分析。

## 流程模块

- `0data_prepare/`：数据下载与预处理（`.py`/`.sh` 脚本 4 个）。
- `1assmbly/`：HiFi/Hi-C 组装流程（`.py`/`.sh` 脚本 4 个）。
- `2SV_identify/`：SV 与倒位识别（`.py`/`.sh` 脚本 12 个）。
- `3Genome_annotation/`：重复序列处理与基因组注释（`.py`/`.sh` 脚本 4 个）。
- `4Post_analysis/`：后续分析与进化统计（`.py`/`.sh` 脚本 13 个）。

## 典型使用顺序

1. 在 `0data_prepare/` 完成数据准备。
2. 在 `1assmbly/` 执行组装流程。
3. 在 `2SV_identify/` 进行倒位识别。
4. 在 `3Genome_annotation/` 完成注释处理。
5. 在 `4Post_analysis/` 进行下游分析。

## 自动化

- GitHub Actions 工作流：`.github/workflows/auto-update-readme.yml`。
- 在脚本变更、手动触发和每周定时任务时自动更新 README。

## 许可证

详见 [LICENSE](LICENSE)。
