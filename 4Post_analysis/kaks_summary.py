#!/usr/bin/env python3
"""
summarize_kaks.py

扫描 ./kaks_re/*.kaks.csv，统计每个样本的 Ka, Ks, Ka/Ks 的 mean, median，
并统计 95th-percentile 以上的 outliers 数量。
过滤规则：
 - 若 Ks == 0 则无法计算 Ka/Ks（计为被排除）
 - 仅对 Ka/Ks <= 2 的记录进行统计（即保留 Ka/Ks<=2）
脚本会输出每个样本的 summary，并生成 kaks_summary.csv
"""

import os
import glob
import re
import pandas as pd
import numpy as np

INPUT_DIR = "kaks_re"
GLOB_PATTERN = os.path.join(INPUT_DIR, "*.kaks.csv")
OUT_CSV = "kaks_summary.csv"

def normalize_colname(s):
    s = s.strip()
    s = re.sub(r'[\s\-]+', '_', s)
    s = re.sub(r'[()/]', '', s)
    return s.lower()

def find_col(df, candidates):
    # candidates: list of possible name fragments, return first matching column name or None
    cols = {normalize_colname(c): c for c in df.columns}
    for cand in candidates:
        nc = normalize_colname(cand)
        # exact or partial matches
        for ncol, orig in cols.items():
            if nc == ncol or nc in ncol or ncol in nc:
                return orig
    return None

summaries = []

files = sorted(glob.glob(GLOB_PATTERN))
if not files:
    raise SystemExit(f"No files found with pattern {GLOB_PATTERN}. Put your .kaks.csv files under '{INPUT_DIR}'.")

for fp in files:
    sample_name = os.path.basename(fp).replace(".kaks.csv", "")
    try:
        # try auto-detect sep (works for comma, tab)
        df = pd.read_csv(fp, sep=None, engine='python')
    except Exception as e:
        print(f"[WARN] Failed to read {fp} with auto sep: {e}. Trying comma fallback.")
        df = pd.read_csv(fp, sep=",", engine='python')

    # normalize column names mapping
    df.columns = [c.strip() for c in df.columns]

    # find relevant columns
    ka_col = find_col(df, ["Ka", "ka", "Nonsynonymous", "dN"])
    ks_col = find_col(df, ["Ks", "ks", "Synonymous", "dS"])
    kak_col = find_col(df, ["Ka/Ks", "Ka_Ks", "ka/ks", "omega", "dN/dS"])

    # read numeric columns robustly
    def to_float_series(col):
        if col is None:
            return pd.Series([np.nan]*len(df))
        s = pd.to_numeric(df[col].replace(['Inf','inf','Infinity','NA','-',''], np.nan), errors='coerce')
        return s

    ka_s = to_float_series(ka_col)
    ks_s = to_float_series(ks_col)
    if kak_col is not None:
        kak_s = to_float_series(kak_col)
    else:
        # compute if possible
        kak_s = pd.Series([np.nan]*len(df))
        with np.errstate(divide='ignore', invalid='ignore'):
            mask = ks_s.notna() & (ks_s != 0)
            kak_s.loc[mask] = ka_s.loc[mask] / ks_s.loc[mask]

    # counts
    total_records = len(df)
    ks_zero_or_na = ((ks_s == 0) | ks_s.isna()).sum()
    # For calculation, treat ks==0 as NaN in kak_s (already is)
    # Now apply the filter: keep only kak_s <= 2 and kak_s not null
    kept_mask = kak_s.notna() & (kak_s <= 2)
    kept_n = kept_mask.sum()
    excluded_by_kak_gt2 = ((kak_s.notna()) & (kak_s > 2)).sum()
    excluded_for_nan_kak = kak_s.isna().sum()  # includes ks==0 or missing ka/ks
    # but excluded_for_nan_kak includes rows where kak missing even if ks>0 but ka missing; that's OK

    # compute stats on kept set
    ka_kept = ka_s[kept_mask]
    ks_kept = ks_s[kept_mask]
    kak_kept = kak_s[kept_mask]

    def summarystats(series):
        if series.dropna().empty:
            return (np.nan, np.nan, np.nan)  # mean, median, 95th threshold
        mean = float(series.mean())
        median = float(series.median())
        p95 = float(np.percentile(series.dropna(), 95))
        return (mean, median, p95)

    ka_mean, ka_median, ka_p95 = summarystats(ka_kept)
    ks_mean, ks_median, ks_p95 = summarystats(ks_kept)
    kak_mean, kak_median, kak_p95 = summarystats(kak_kept)

    # outlier counts = number strictly > 95th percentile (computed on kept data)
    ka_outliers = int((ka_kept > ka_p95).sum()) if not np.isnan(ka_p95) else 0
    ks_outliers = int((ks_kept > ks_p95).sum()) if not np.isnan(ks_p95) else 0
    kak_outliers = int((kak_kept > kak_p95).sum()) if not np.isnan(kak_p95) else 0

    summary = {
        "sample": sample_name,
        "total_records": int(total_records),
        "kept_records_ka_ks_ratio_le_2": int(kept_n),
        "excluded_by_kak_gt_2": int(excluded_by_kak_gt2),
        "excluded_kak_missing_or_ks_zero": int(excluded_for_nan_kak),
        "ks_zero_or_missing_count": int(ks_zero_or_na),

        "ka_mean_kept": ka_mean,
        "ka_median_kept": ka_median,
        "ka_95pct_threshold_kept": ka_p95,
        "ka_outliers_gt_95pct": ka_outliers,

        "ks_mean_kept": ks_mean,
        "ks_median_kept": ks_median,
        "ks_95pct_threshold_kept": ks_p95,
        "ks_outliers_gt_95pct": ks_outliers,

        "kak_mean_kept": kak_mean,
        "kak_median_kept": kak_median,
        "kak_95pct_threshold_kept": kak_p95,
        "kak_outliers_gt_95pct": kak_outliers,
    }

    summaries.append(summary)
    print(f"[OK] {sample_name}  total={total_records}  kept={kept_n}  excl(>2)={excluded_by_kak_gt2}  ks==0_or_missing={ks_zero_or_na}")

# save summary table
out_df = pd.DataFrame(summaries)
out_df.to_csv(OUT_CSV, index=False)
print(f"\nWrote summary for {len(summaries)} samples to {OUT_CSV}")