#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute simple inversion-vs-background Ka/Ks shift metrics.

Inputs
------
1) inv_gene_table: tab/csv with at least
      sample_id, inv_id, gene_id
2) kaks_dir: directory containing one genome-wide Ka/Ks table per sample

Expected Ka/Ks table format
---------------------------
The script auto-detects delimiter and common column names.
It expects one gene/transcript identifier column and one Ka/Ks (omega) column.
Examples of supported columns:
  Gene, gene_id, GeneID, transcript_id
  Ka/Ks, omega, kaks, KaKs

Key behavior
------------
- Transcript IDs such as g756.t1 / g756.t2 are collapsed to gene IDs (g756)
- Multiple transcript-level Ka/Ks rows per gene are collapsed to a single
  gene-level value using the median
- Background is defined as all genes in the sample Ka/Ks table excluding
  all inversion genes from that sample
- For each inversion, compute:
    * inversion median Ka/Ks
    * inversion p90 Ka/Ks
    * background median Ka/Ks
    * background p90 Ka/Ks
    * kaks_shift_ratio_median = inv_median / bg_median
    * kaks_shift_ratio_p90    = inv_p90 / bg_p90
    * corresponding delta values

This is intentionally simple and direct for downstream figure comparisons.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


GENE_COL_CANDIDATES = [
    "gene_id", "gene", "Gene", "GeneID", "geneID", "transcript_id",
    "transcript", "Sequence", "seq_id", "ID"
]
OMEGA_COL_CANDIDATES = [
    "Ka/Ks", "omega", "Omega", "kaks", "KaKs", "ka_ks", "dN/dS"
]


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def sniff_delimiter(path: Path, default: str = "\t") -> str:
    try:
        with path.open("r", encoding="utf-8", errors="ignore") as fh:
            sample = fh.read(4096)
        dialect = csv.Sniffer().sniff(sample, delimiters=[",", "\t", ";"])
        return dialect.delimiter
    except Exception:
        with path.open("r", encoding="utf-8", errors="ignore") as fh:
            first = fh.readline()
        if first.count("\t") >= first.count(",") and first.count("\t") > 0:
            return "\t"
        if first.count(",") > 0:
            return ","
        return default


_TRANSCRIPT_PATTERNS = [
    re.compile(r"(\.t\d+)$", re.IGNORECASE),
    re.compile(r"(\.mrna\d+)$", re.IGNORECASE),
    re.compile(r"(-T\d+)$", re.IGNORECASE),
    re.compile(r"(_T\d+)$", re.IGNORECASE),
    re.compile(r"(\.isoform\d+)$", re.IGNORECASE),
]


def strip_transcript_suffix(gene_id: str) -> str:
    s = str(gene_id).strip()
    for pat in _TRANSCRIPT_PATTERNS:
        s = pat.sub("", s)
    return s


def choose_existing_column(df: pd.DataFrame, candidates: Sequence[str]) -> Optional[str]:
    cols = set(df.columns)
    for c in candidates:
        if c in cols:
            return c
    lower_to_orig = {str(c).lower(): c for c in df.columns}
    for c in candidates:
        key = c.lower()
        if key in lower_to_orig:
            return lower_to_orig[key]
    return None


def read_table_auto(path: Path) -> pd.DataFrame:
    delim = sniff_delimiter(path)
    try:
        return pd.read_csv(path, sep=delim, dtype=str, low_memory=False)
    except Exception:
        return pd.read_csv(path, sep=delim, dtype=str, engine="python")


def parse_numeric_series(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s.astype(str).str.strip(), errors="coerce")


def normalize_sample_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "", str(s)).lower()


def build_sample_file_map(kaks_dir: Path, suffix: Optional[str]) -> Dict[str, Path]:
    files: List[Path] = []
    if suffix:
        files = sorted(kaks_dir.glob(f"*{suffix}"))
    else:
        files = sorted([p for p in kaks_dir.iterdir() if p.is_file()])

    fmap: Dict[str, Path] = {}
    for fp in files:
        stem = fp.name
        if suffix and stem.endswith(suffix):
            stem = stem[: -len(suffix)]
        else:
            stem = fp.stem
        fmap[normalize_sample_name(stem)] = fp
    return fmap


def load_kaks_gene_table(
    path: Path,
    omega_min: float = 0.0,
    omega_max: float = 10.0,
    collapse_fun: str = "median",
) -> pd.DataFrame:
    raw = read_table_auto(path)
    gene_col = choose_existing_column(raw, GENE_COL_CANDIDATES)
    omega_col = choose_existing_column(raw, OMEGA_COL_CANDIDATES)

    if gene_col is None or omega_col is None:
        raise ValueError(
            f"{path.name}: cannot find gene and/or Ka/Ks column. "
            f"Columns present: {list(raw.columns)}"
        )

    df = raw[[gene_col, omega_col]].copy()
    df.columns = ["gene_raw", "omega"]
    df["gene_id"] = df["gene_raw"].map(strip_transcript_suffix)
    df["omega"] = parse_numeric_series(df["omega"])
    df = df.dropna(subset=["gene_id", "omega"])
    df = df[np.isfinite(df["omega"])]
    df = df[df["omega"] >= omega_min]
    if omega_max is not None and omega_max > 0:
        df = df[df["omega"] <= omega_max]

    if collapse_fun == "mean":
        agg = df.groupby("gene_id", as_index=False)["omega"].mean()
    else:
        agg = df.groupby("gene_id", as_index=False)["omega"].median()

    return agg


def safe_quantile(x: np.ndarray, q: float) -> float:
    if x.size == 0:
        return float("nan")
    return float(np.quantile(x, q, method="linear" if hasattr(np, "quantile") else "midpoint"))


def safe_ratio(num: float, den: float) -> float:
    if pd.isna(num) or pd.isna(den) or den == 0:
        return float("nan")
    return float(num / den)


def summarize_inversion(inv_values: np.ndarray, bg_values: np.ndarray) -> Dict[str, float]:
    inv_median = float(np.median(inv_values)) if inv_values.size else float("nan")
    inv_p90 = safe_quantile(inv_values, 0.90) if inv_values.size else float("nan")
    bg_median = float(np.median(bg_values)) if bg_values.size else float("nan")
    bg_p90 = safe_quantile(bg_values, 0.90) if bg_values.size else float("nan")

    return {
        "inv_kaks_median": inv_median,
        "inv_kaks_p90": inv_p90,
        "bg_kaks_median": bg_median,
        "bg_kaks_p90": bg_p90,
        "kaks_shift_ratio_median": safe_ratio(inv_median, bg_median),
        "kaks_shift_ratio_p90": safe_ratio(inv_p90, bg_p90),
        "kaks_shift_delta_median": (inv_median - bg_median) if np.isfinite(inv_median) and np.isfinite(bg_median) else float("nan"),
        "kaks_shift_delta_p90": (inv_p90 - bg_p90) if np.isfinite(inv_p90) and np.isfinite(bg_p90) else float("nan"),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description="Simple inversion-vs-background Ka/Ks shift metrics")
    ap.add_argument("--inv-genes", required=True, help="TSV/CSV with sample_id, inv_id, gene_id")
    ap.add_argument("--kaks-dir", required=True, help="Directory of per-sample genome-wide Ka/Ks tables")
    ap.add_argument("--output", required=True, help="Output TSV path")
    ap.add_argument("--suffix", default=None, help="Optional suffix for Ka/Ks files, e.g. .csv")
    ap.add_argument("--omega-min", type=float, default=0.0, help="Minimum allowed Ka/Ks (default: 0)")
    ap.add_argument("--omega-max", type=float, default=10.0, help="Maximum allowed Ka/Ks; set <=0 to disable")
    ap.add_argument("--collapse-fun", choices=["median", "mean"], default="median", help="Collapse transcript-level Ka/Ks to gene-level")
    args = ap.parse_args()

    inv_path = Path(args.inv_genes)
    kaks_dir = Path(args.kaks_dir)
    out_path = Path(args.output)

    inv_df = read_table_auto(inv_path)
    required = {"sample_id", "inv_id", "gene_id"}
    if not required.issubset(inv_df.columns):
        raise ValueError(f"inv_gene_table must contain columns {sorted(required)}; got {list(inv_df.columns)}")

    inv_df = inv_df[["sample_id", "inv_id", "gene_id"]].copy()
    inv_df["sample_id"] = inv_df["sample_id"].astype(str).str.strip()
    inv_df["inv_id"] = inv_df["inv_id"].astype(str).str.strip()
    inv_df["gene_id"] = inv_df["gene_id"].astype(str).map(strip_transcript_suffix)
    inv_df = inv_df.dropna().drop_duplicates()

    sample_file_map = build_sample_file_map(kaks_dir, args.suffix)
    omega_max = None if args.omega_max is None or args.omega_max <= 0 else args.omega_max

    results: List[Dict[str, object]] = []

    for sample_id, sub in inv_df.groupby("sample_id", sort=True):
        norm_sample = normalize_sample_name(sample_id)
        fp = sample_file_map.get(norm_sample)
        if fp is None:
            # try contains-based fallback
            candidates = [p for key, p in sample_file_map.items() if norm_sample in key or key in norm_sample]
            fp = candidates[0] if len(candidates) == 1 else None

        if fp is None:
            eprint(f"[WARN] No Ka/Ks file matched sample_id={sample_id}")
            for inv_id, inv_sub in sub.groupby("inv_id"):
                results.append({
                    "sample_id": sample_id,
                    "inv_id": inv_id,
                    "status": "kaks_file_not_found",
                    "n_genes_input": int(inv_sub["gene_id"].nunique()),
                })
            continue

        try:
            gene_tbl = load_kaks_gene_table(
                fp,
                omega_min=args.omega_min,
                omega_max=omega_max,
                collapse_fun=args.collapse_fun,
            )
        except Exception as exc:
            eprint(f"[WARN] Failed to load {fp.name}: {exc}")
            for inv_id, inv_sub in sub.groupby("inv_id"):
                results.append({
                    "sample_id": sample_id,
                    "inv_id": inv_id,
                    "status": f"kaks_load_failed:{type(exc).__name__}",
                    "n_genes_input": int(inv_sub["gene_id"].nunique()),
                })
            continue

        gene_to_omega = dict(zip(gene_tbl["gene_id"], gene_tbl["omega"]))
        sample_all_inv_genes = set(sub["gene_id"].unique())
        bg_tbl = gene_tbl[~gene_tbl["gene_id"].isin(sample_all_inv_genes)].copy()
        bg_values = bg_tbl["omega"].to_numpy(dtype=float)

        if bg_values.size == 0:
            eprint(f"[WARN] Empty background after excluding inversion genes for sample_id={sample_id}")
            for inv_id, inv_sub in sub.groupby("inv_id"):
                results.append({
                    "sample_id": sample_id,
                    "inv_id": inv_id,
                    "status": "empty_background",
                    "n_genes_input": int(inv_sub["gene_id"].nunique()),
                })
            continue

        bg_median = float(np.median(bg_values))
        bg_p90 = safe_quantile(bg_values, 0.90)
        bg_n = int(bg_values.size)
        total_gene_n = int(gene_tbl.shape[0])

        for inv_id, inv_sub in sub.groupby("inv_id", sort=False):
            inv_gene_set = list(dict.fromkeys(inv_sub["gene_id"].tolist()))
            found = [g for g in inv_gene_set if g in gene_to_omega]
            missing = [g for g in inv_gene_set if g not in gene_to_omega]
            inv_values = np.array([gene_to_omega[g] for g in found], dtype=float)

            row: Dict[str, object] = {
                "sample_id": sample_id,
                "inv_id": inv_id,
                "status": "ok" if inv_values.size > 0 else "no_genes_found_in_kaks",
                "kaks_file": fp.name,
                "n_genes_input": len(inv_gene_set),
                "n_genes_found": int(inv_values.size),
                "n_genes_missing": len(missing),
                "sample_total_genes_in_kaks": total_gene_n,
                "sample_background_genes": bg_n,
                "bg_kaks_median": bg_median,
                "bg_kaks_p90": bg_p90,
            }

            if inv_values.size > 0:
                row.update(summarize_inversion(inv_values, bg_values))
            else:
                row.update({
                    "inv_kaks_median": np.nan,
                    "inv_kaks_p90": np.nan,
                    "kaks_shift_ratio_median": np.nan,
                    "kaks_shift_ratio_p90": np.nan,
                    "kaks_shift_delta_median": np.nan,
                    "kaks_shift_delta_p90": np.nan,
                })
            results.append(row)

    out_df = pd.DataFrame(results)

    # Column order
    wanted = [
        "sample_id", "inv_id", "status", "kaks_file",
        "n_genes_input", "n_genes_found", "n_genes_missing",
        "sample_total_genes_in_kaks", "sample_background_genes",
        "bg_kaks_median", "bg_kaks_p90",
        "inv_kaks_median", "inv_kaks_p90",
        "kaks_shift_ratio_median", "kaks_shift_ratio_p90",
        "kaks_shift_delta_median", "kaks_shift_delta_p90",
    ]
    existing = [c for c in wanted if c in out_df.columns]
    remainder = [c for c in out_df.columns if c not in existing]
    out_df = out_df[existing + remainder]
    out_df.to_csv(out_path, sep="\t", index=False)
    eprint(f"[INFO] Wrote {out_df.shape[0]} rows to {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
