#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import math
from pathlib import Path
from typing import Optional

import pandas as pd


def detect_sep(path: str) -> str:
    """
    Detect delimiter for KaKs table.
    Tries tab first (most likely), then comma, semicolon, whitespace.
    """
    candidates = ["\t", ",", ";"]
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        # sample some non-empty lines
        lines = []
        for _ in range(30):
            line = f.readline()
            if not line:
                break
            line = line.strip("\n")
            if line.strip():
                lines.append(line)
        if not lines:
            return "\t"

    best_sep = "\t"
    best_score = -1
    for sep in candidates:
        counts = [len(l.split(sep)) for l in lines]
        avg = sum(counts) / len(counts)
        if avg > best_score:
            best_score = avg
            best_sep = sep

    # If still looks like 1 column, fall back to any whitespace
    if best_score <= 1.2:
        return r"\s+"
    return best_sep


def find_col(df: pd.DataFrame, target: str) -> Optional[str]:
    """
    Find a column name in df that matches target case-insensitively after stripping.
    """
    target_lc = target.strip().lower()
    for c in df.columns:
        if c.strip().lower() == target_lc:
            return c
    return None


def main():
    ap = argparse.ArgumentParser(
        description="Hard-filter KaKs table by Ks range and alignment length."
    )
    ap.add_argument("--infile", required=True, help="Input KaKs table (TSV/CSV).")
    ap.add_argument("--outfile", required=True, help="Output filtered table (TSV).")
    ap.add_argument("--ks_min", type=float, default=0.01, help="Minimum Ks (default 0.01).")
    ap.add_argument("--ks_max", type=float, default=2.0, help="Maximum Ks (default 2).")
    ap.add_argument("--len_min", type=int, default=300, help="Minimum Length in bp (default 300).")
    ap.add_argument("--out_sep", default="\t", help="Output delimiter: default TAB.")
    args = ap.parse_args()

    infile = args.infile
    outfile = args.outfile

    sep = detect_sep(infile)
    df = pd.read_csv(infile, sep=sep, engine="python", dtype=str, comment=None)
    # strip column whitespace
    df.columns = [c.strip() for c in df.columns]

    # locate required columns
    col_gene = find_col(df, "Gene")
    col_ka   = find_col(df, "Ka")
    col_ks   = find_col(df, "Ks")
    col_len  = find_col(df, "Length")

    missing = [name for name, col in [("Gene", col_gene), ("Ka", col_ka), ("Ks", col_ks), ("Length", col_len)] if col is None]
    if missing:
        raise SystemExit(f"[ERROR] Missing required columns: {', '.join(missing)}. Found columns: {list(df.columns)}")

    # coerce numeric for Ka/Ks/Length
    df[col_ka]  = pd.to_numeric(df[col_ka], errors="coerce")
    df[col_ks]  = pd.to_numeric(df[col_ks], errors="coerce")
    df[col_len] = pd.to_numeric(df[col_len], errors="coerce")

    n0 = df.shape[0]

    # drop NA in required numeric fields
    m = df[col_ka].notna() & df[col_ks].notna() & df[col_len].notna()
    df1 = df.loc[m].copy()
    n1 = df1.shape[0]

    # finite check (handles inf/-inf)
    m = df1[col_ka].map(math.isfinite) & df1[col_ks].map(math.isfinite) & df1[col_len].map(math.isfinite)
    df2 = df1.loc[m].copy()
    n2 = df2.shape[0]

    # length filter
    m = df2[col_len] >= args.len_min
    df3 = df2.loc[m].copy()
    n3 = df3.shape[0]

    # Ks range filter + Ks > 0 (implicit stability)
    m = (df3[col_ks] >= args.ks_min) & (df3[col_ks] <= args.ks_max)
    df4 = df3.loc[m].copy()
    n4 = df4.shape[0]

    # write output
    out_path = Path(outfile)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df4.to_csv(outfile, sep=args.out_sep, index=False)

    # report filtering stats to stderr
    print(f"[INFO] Input rows: {n0}", flush=True)
    print(f"[INFO] After drop NA (Ka/Ks/Length): {n1}", flush=True)
    print(f"[INFO] After finite check: {n2}", flush=True)
    print(f"[INFO] After Length >= {args.len_min}: {n3}", flush=True)
    print(f"[INFO] After {args.ks_min} <= Ks <= {args.ks_max}: {n4}", flush=True)
    print(f"[INFO] Wrote filtered table: {outfile}", flush=True)


if __name__ == "__main__":
    main()