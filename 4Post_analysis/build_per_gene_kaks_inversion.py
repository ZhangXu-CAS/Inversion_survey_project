#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import math
import os
import re
import sys
from typing import Dict, Optional, Tuple, Iterable


def safe_float(s: str) -> Optional[float]:
    if s is None:
        return None
    s = str(s).strip()
    if s == "" or s.upper() in {"NA", "N/A", "NAN", "NONE", "NULL"}:
        return None
    try:
        x = float(s)
        if math.isinf(x) or math.isnan(x) or x < 0:
            return None
        return x
    except Exception:
        return None


def normalize(s: str) -> str:
    return "".join(ch.lower() for ch in s if ch.isalnum())


def detect_column(fieldnames, candidates: Tuple[str, ...]) -> Optional[str]:
    if not fieldnames:
        return None
    norm_map = {normalize(k): k for k in fieldnames}
    for c in candidates:
        key = normalize(c)
        if key in norm_map:
            return norm_map[key]
    return None


def sniff_delimiter(path: str, default: str = ",") -> str:
    with open(path, "r", encoding="utf-8", newline="") as f:
        head = f.read(4096)
    first = head.splitlines()[0] if head else ""
    if "\t" in first:
        return "\t"
    try:
        dialect = csv.Sniffer().sniff(head, delimiters=[",", "\t", ";", "|"])
        return dialect.delimiter
    except Exception:
        return default


def normalize_gene_id(g: str, mode: str = "gene") -> str:
    """
    mode="gene":  g756.t1 -> g756
    mode="full":  keep full id (g756.t1)
    Also strips common prefixes/spaces.
    """
    if g is None:
        return ""
    g = str(g).strip()

    # strip prefixes like "gene:" "ID=" etc if exist
    g = re.sub(r"^(gene:|Gene:|ID=|gene_id=)", "", g).strip()

    if mode == "full":
        return g

    # mode == "gene": take before first dot, e.g. g1000.t1 -> g1000
    # also handle cases like g1000.t1;something (shouldn't happen but safe)
    g = g.split(";")[0].strip()
    base = g.split(".", 1)[0].strip()
    return base


def load_gene_omega_map(
    gene_kaks_path: str,
    gene_col_candidates: Tuple[str, ...],
    omega_col_candidates: Tuple[str, ...],
    gene_id_mode: str = "gene",
) -> Dict[str, Optional[float]]:
    """
    Build map: normalized_gene_id -> omega
    If multiple transcripts map to same gene, keep the FIRST non-NA omega encountered.
    """
    delim = sniff_delimiter(gene_kaks_path, default=",")

    gene2omega: Dict[str, Optional[float]] = {}

    with open(gene_kaks_path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=delim)
        if reader.fieldnames is None:
            raise ValueError(f"No header found in {gene_kaks_path}")

        gene_col = detect_column(reader.fieldnames, gene_col_candidates)
        omega_col = detect_column(reader.fieldnames, omega_col_candidates)

        if gene_col is None:
            raise ValueError(
                f"Cannot find gene column in {gene_kaks_path}. "
                f"Delimiter detected: {repr(delim)}. "
                f"Tried: {gene_col_candidates}. Available: {reader.fieldnames}"
            )
        if omega_col is None:
            raise ValueError(
                f"Cannot find omega/KaKs column in {gene_kaks_path}. "
                f"Delimiter detected: {repr(delim)}. "
                f"Tried: {omega_col_candidates}. Available: {reader.fieldnames}"
            )

        for row in reader:
            raw_gene = row.get(gene_col)
            gene = normalize_gene_id(raw_gene, mode=gene_id_mode)
            if not gene:
                continue

            omega = safe_float(row.get(omega_col))

            if gene not in gene2omega:
                gene2omega[gene] = omega
            else:
                # upgrade None -> valid
                if gene2omega[gene] is None and omega is not None:
                    gene2omega[gene] = omega

    return gene2omega


def iter_inversion_gene_pairs(
    inv_tsv: str,
    sample_col: str,
    inv_col: str,
    genes_col: str,
    genes_sep: str = ";",
    inv_gene_id_mode: str = "gene",
) -> Iterable[Tuple[str, str, str, str]]:
    """
    Yield (sample_id, inv_id, raw_gene_id, normalized_gene_id)
    """
    with open(inv_tsv, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"No header found in {inv_tsv}")

        for c in (sample_col, inv_col, genes_col):
            if c not in reader.fieldnames:
                raise ValueError(
                    f"Missing required column '{c}' in {inv_tsv}. "
                    f"Available columns: {reader.fieldnames}"
                )

        for row in reader:
            sample_id = (row.get(sample_col) or "").strip()
            inv_id = (row.get(inv_col) or "").strip()
            genes_raw = (row.get(genes_col) or "").strip()
            if not sample_id or not inv_id or not genes_raw:
                continue

            for g in genes_raw.split(genes_sep):
                raw_gene = g.strip()
                if not raw_gene:
                    continue
                norm_gene = normalize_gene_id(raw_gene, mode=inv_gene_id_mode)
                yield sample_id, inv_id, raw_gene, norm_gene


def main():
    ap = argparse.ArgumentParser(
        description="Build per_gene_kaks_inversion.csv from inv_genes_per_inversion.tsv and per-sample gene Ka/Ks files; handles gXXX.t1 vs gXXX."
    )
    ap.add_argument("--inv_tsv", required=True, help="inv_genes_per_inversion.tsv (TSV)")
    ap.add_argument("--kaks_dir", required=True, help="Directory containing <sample>.gene_kaks.csv (often TSV)")
    ap.add_argument("--out_csv", required=True, help="Output CSV path")

    ap.add_argument("--sample_col", default="sample")
    ap.add_argument("--inv_col", default="inv_id")
    ap.add_argument("--genes_col", default="genes_inside_ids")
    ap.add_argument("--genes_sep", default=";")

    ap.add_argument("--gene_col_candidates", default="Gene,gene,gene_id,GeneID,ID")
    ap.add_argument("--omega_col_candidates", default="Ka/Ks,KaKs,omega,Ka_Ks,KaKsRatio,KaKs_ratio")

    # gene ID normalization controls
    ap.add_argument("--kaks_gene_id_mode", default="gene", choices=["gene", "full"],
                    help="How to normalize gene IDs from gene_kaks file: gene (g1.t1->g1) or full (keep). Default: gene.")
    ap.add_argument("--inv_gene_id_mode", default="gene", choices=["gene", "full"],
                    help="How to normalize gene IDs from inversion TSV: gene (g1.t1->g1) or full. Default: gene.")

    ap.add_argument("--missing_kaks", default="warn", choices=["warn", "skip", "error"],
                    help="What to do if a sample's gene_kaks file is missing.")

    # diagnostics
    ap.add_argument("--report_top_missing", type=int, default=20,
                    help="Report top N missing gene IDs per sample (default 20).")
    args = ap.parse_args()

    gene_col_candidates = tuple(x.strip() for x in args.gene_col_candidates.split(",") if x.strip())
    omega_col_candidates = tuple(x.strip() for x in args.omega_col_candidates.split(",") if x.strip())

    cache: Dict[str, Dict[str, Optional[float]]] = {}

    n_rows = 0
    n_missing_sample_file = 0
    n_missing_gene = 0

    # per-sample missing gene counters (for debug)
    missing_gene_counter: Dict[Tuple[str, str], int] = {}

    os.makedirs(os.path.dirname(os.path.abspath(args.out_csv)) or ".", exist_ok=True)
    with open(args.out_csv, "w", newline="", encoding="utf-8") as out_f:
        w = csv.writer(out_f)
        w.writerow(["sample_id", "inv_id", "gene_id", "omega"])

        for sample_id, inv_id, raw_gene, norm_gene in iter_inversion_gene_pairs(
            args.inv_tsv, args.sample_col, args.inv_col, args.genes_col,
            genes_sep=args.genes_sep, inv_gene_id_mode=args.inv_gene_id_mode
        ):
            if sample_id not in cache:
                gene_kaks_path = os.path.join(args.kaks_dir, f"{sample_id}.gene_kaks.csv")
                if not os.path.exists(gene_kaks_path):
                    n_missing_sample_file += 1
                    msg = f"[missing] {gene_kaks_path}"
                    if args.missing_kaks == "error":
                        raise FileNotFoundError(msg)
                    elif args.missing_kaks == "warn":
                        print(msg, file=sys.stderr)
                    cache[sample_id] = {}
                else:
                    cache[sample_id] = load_gene_omega_map(
                        gene_kaks_path,
                        gene_col_candidates,
                        omega_col_candidates,
                        gene_id_mode=args.kaks_gene_id_mode
                    )

            omega = cache[sample_id].get(norm_gene)
            if omega is None:
                n_missing_gene += 1
                missing_gene_counter[(sample_id, norm_gene)] = missing_gene_counter.get((sample_id, norm_gene), 0) + 1
                w.writerow([sample_id, inv_id, norm_gene, "NA"])
            else:
                w.writerow([sample_id, inv_id, norm_gene, f"{omega:.10g}"])
            n_rows += 1

    print(f"[done] wrote {n_rows} rows -> {args.out_csv}", file=sys.stderr)
    print(f"[stats] missing sample gene_kaks files: {n_missing_sample_file}", file=sys.stderr)
    print(f"[stats] missing gene omega entries: {n_missing_gene}", file=sys.stderr)

    # report top missing genes (helpful for debugging mismatched naming)
    if args.report_top_missing > 0 and missing_gene_counter:
        # aggregate by sample then list most frequent missing genes
        from collections import defaultdict

        by_sample = defaultdict(list)
        for (sid, gid), cnt in missing_gene_counter.items():
            by_sample[sid].append((cnt, gid))
        for sid, lst in by_sample.items():
            lst.sort(reverse=True)
            top = lst[: args.report_top_missing]
            print(f"[debug] sample={sid} top_missing_genes (count, gene): {top}", file=sys.stderr)


if __name__ == "__main__":
    main()