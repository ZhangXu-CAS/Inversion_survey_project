#!/bin/bash
#SBATCH --job-name=KaKs_1031
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=80G
#SBATCH --time=6:00:00
#SBATCH --account=def-rieseber
#SBATCH --array=0-101
#SBATCH --output=log/KaKs_%A_%a.out

set -euo pipefail
# --- config: change if needed ---
PROJECT_ROOT=/project/6004758/zhangxu/Gene_Ks
RESULTS_ROOT=${PROJECT_ROOT}   # where important outputs will be copied back
MERGE_AXT=${PROJECT_ROOT}/merge_codonfas_to_axt.py
PY_HELPER=${PROJECT_ROOT}/backtranslate_codon_align.py
PAV_CAL=${PROJECT_ROOT}/pav_from_blast.py

# number of parallel alignment jobs for step 4 (on-node parallelism)
PAR_JOBS=20  # set via env if you want: export PAR_JOBS=16

# pick thread count consistent with SLURM
THREADS=20

# Tools
source /home/zhangxu/miniconda3/bin/activate kaks
module load muscle/5.1.0
module load blast+/2.14.1

KAKS_CMD=$(command -v KaKs_Calculator || true)
MUSCLE_CMD=$(command -v muscle || true)
BLASTP_CMD=$(command -v blastp || true)
MAKEBLASTDB_CMD=$(command -v makeblastdb || true)

# --- sample id from array ---
cd "$PROJECT_ROOT"
id=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" ks_ids2)
echo "[$(date)] START sample: $id (SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID)"

# input files on project storage
HAP1_PROT=/project/6004758/zhangxu/TE_anno/${id}/braker_out1/braker.aa
HAP2_PROT=/project/6004758/zhangxu/TE_anno/${id}/braker_out2/braker.aa
HAP1_CDS=/project/6004758/zhangxu/TE_anno/${id}/braker_out1/braker.cds
HAP2_CDS=/project/6004758/zhangxu/TE_anno/${id}/braker_out2/braker.cds


# node-local tmp dir on host (created here)
TMP_BASE="${SLURM_TMPDIR:-/tmp}"
WORKDIR=$(mktemp -d "${TMP_BASE}/kaks_work_${id}.XXXXX")
echo "WORKDIR (host) = $WORKDIR"
cd $WORKDIR
# cleanup and copy-back handler (will be run on exit)
function copy_back_and_cleanup {
    echo "[$(date)] trap: copying back results (if any) and cleaning up..."
    mkdir -p "${RESULTS_ROOT}/${id}"
    # copy important results if present
    cp -a --preserve=mode,timestamps "${WORKDIR}/hap1_vs_hap2.blast" "${RESULTS_ROOT}/${id}/" 2>/dev/null || true
    cp -a --preserve=mode,timestamps "${WORKDIR}/hap2_vs_hap1.blast" "${RESULTS_ROOT}/${id}/" 2>/dev/null || true
    if [[ -d "${WORKDIR}/pav_blast_out" ]]; then
        rsync -a "${WORKDIR}/pav_blast_out" "${RESULTS_ROOT}/${id}/" || cp -a "${WORKDIR}/pav_blast_out" "${RESULTS_ROOT}/${id}/"
    fi
    if [[ -s "${WORKDIR}/${id}.merged.axt" ]]; then
        cp -a "${WORKDIR}/${id}.merged.axt" "${RESULTS_ROOT}/${id}/" 2>/dev/null || true
    fi
    if [[ -s "${WORKDIR}/${id}.kaks.csv" ]]; then
        cp -a "${WORKDIR}/${id}.kaks.csv" "${RESULTS_ROOT}/${id}/" 2>/dev/null || true
    fi
    if [[ -s "${WORKDIR}/${id}.intermediates.tar.gz" ]]; then
        cp -a "${WORKDIR}/${id}.intermediates.tar.gz" "${RESULTS_ROOT}/${id}/" 2>/dev/null || true
    fi
    
}
trap copy_back_and_cleanup EXIT


# ---- copy inputs to local node (avoid NFS heavy IO during many small-file ops) ----
echo "Copying input fasta files to node-local disk..."
cp -a "$HAP1_PROT" "$WORKDIR/hap1.prot.fa" || { echo "Missing $HAP1_PROT"; exit 1; }
cp -a "$HAP2_PROT" "$WORKDIR/hap2.prot.fa" || { echo "Missing $HAP2_PROT"; exit 1; }
cp -a "$HAP1_CDS" "$WORKDIR/hap1.cds.fa" || { echo "Missing $HAP1_CDS"; exit 1; }
cp -a "$HAP2_CDS" "$WORKDIR/hap2.cds.fa" || { echo "Missing $HAP2_CDS"; exit 1; }

# rename local copies to simple names
LOCAL_HAP1_PROT="$WORKDIR/hap1.prot.fa"
LOCAL_HAP2_PROT="$WORKDIR/hap2.prot.fa"
LOCAL_HAP1_CDS="$WORKDIR/hap1.cds.fa"
LOCAL_HAP2_CDS="$WORKDIR/hap2.cds.fa"

# ========== 1) makeblastdb & blastp on node local ==========
echo "[1] makeblastdb..."
$MAKEBLASTDB_CMD -dbtype prot -in "$LOCAL_HAP1_PROT" -out hap1_db
$MAKEBLASTDB_CMD -dbtype prot -in "$LOCAL_HAP2_PROT" -out hap2_db

echo "[1] blastp hap1 -> hap2..."
$BLASTP_CMD -query "$LOCAL_HAP1_PROT" -db hap2_db -out hap1_vs_hap2.blast -evalue 1e-5 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -num_threads $THREADS

echo "[1] blastp hap2 -> hap1..."
$BLASTP_CMD -query "$LOCAL_HAP2_PROT" -db hap1_db -out hap2_vs_hap1.blast -evalue 1e-5 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -num_threads $THREADS

# ========== 2) RBH detection (top-hit) on node local ==========
echo "[2] parsing RBH (reciprocal best hits)..."
sort -k1,1 -k12,12nr hap1_vs_hap2.blast | awk '!seen[$1]++{print $1"\t"$2}' > hap1.top.tsv
sort -k1,1 -k12,12nr hap2_vs_hap1.blast | awk '!seen[$1]++{print $1"\t"$2}' > hap2.top.tsv
awk 'NR==FNR {a[$1]=$2; next} { if (a[$2] == $1) print $2"\t"$1 }' hap1.top.tsv hap2.top.tsv > rbh_pairs.tsv
RBH_COUNT=$(wc -l < rbh_pairs.tsv || true)
echo "RBH pairs found: $RBH_COUNT"

# ========== 3) PAV via pav_from_blast.py (local) ==========
echo "[3] PAV stats..."
mkdir -p pav_blast_out
python3 "$PAV_CAL" --hap1-prot "$LOCAL_HAP1_PROT" --hap2-prot "$LOCAL_HAP2_PROT" \
    --b1 hap1_vs_hap2.blast --b2 hap2_vs_hap1.blast --outdir pav_blast_out

# ========== 4) per-pair protein align -> codon align (local) ==========
echo "[4] producing codon alignments for RBH pairs (node-local). Parallel jobs: $PAR_JOBS"
mkdir -p codon_aligns protein_aligns cds_pairs

# function to process one pair (safe)
process_pair() {
    q=$1; s=$2
    # extract proteins
    awk -v q="$q" 'BEGIN{RS=">"; ORS=""} NR>1{split($0,l,"\n"); if (l[1]==q) print ">"$0}' "$LOCAL_HAP1_PROT" > "protein_aligns/${q}_${s}.hap1.prot.fa"
    awk -v q="$s" 'BEGIN{RS=">"; ORS=""} NR>1{split($0,l,"\n"); if (l[1]==q) print ">"$0}' "$LOCAL_HAP2_PROT" > "protein_aligns/${q}_${s}.hap2.prot.fa"
    cat "protein_aligns/${q}_${s}.hap1.prot.fa" "protein_aligns/${q}_${s}.hap2.prot.fa" > "protein_aligns/${q}_${s}.prot.fa"
    # align
    if [[ -n "$MUSCLE_CMD" ]]; then
        $MUSCLE_CMD -align "protein_aligns/${q}_${s}.prot.fa" -output "protein_aligns/${q}_${s}.prot.aln.fa" 
        python3 "$PY_HELPER" "protein_aligns/${q}_${s}.prot.aln.fa" "$LOCAL_HAP1_CDS" "$LOCAL_HAP2_CDS" "$q" "$s" "codon_aligns/${q}_${s}.codon.aln.fa" || echo "backtranslate fail ${q}_${s}" >&2
    else
        echo "muscle missing; skipping ${q}_${s}" >&2
    fi
}

export -f process_pair
export LOCAL_HAP1_PROT LOCAL_HAP2_PROT LOCAL_HAP1_CDS LOCAL_HAP2_CDS PY_HELPER MUSCLE_CMD

# run pairs in parallel if xargs available
awk '{print $1"\t"$2}' hap1.top.tsv | xargs -P "$PAR_JOBS" -n2 bash -c 'process_pair "$1" "$2"' _ 

# ========== 5) merge codon align -> merged.axt and run KaKs (local) ==========
echo "[5] Ka/Ks calculation"
python3 "$MERGE_AXT" codon_aligns "${id}.merged.axt"

if [[ -s "${id}.merged.axt" && -n "$KAKS_CMD" ]]; then
    "$KAKS_CMD" -i "${id}.merged.axt" -o "${id}.kaks.csv" -m "YN" -c 1 || echo "KaKs failed for $id" >&2
else
    echo "No merged AXT or KaKs_Calculator missing; skipping KaKs for $id" >&2
fi

# ========== 6) tar intermediate dirs (single tar) and remove local dirs ==========
echo "[6] packaging intermediates into single tar to reduce file count on shared FS..."
INTERMEDIATE_TAR="${id}.intermediates.tar.gz"
tar -czf "$INTERMEDIATE_TAR" protein_aligns codon_aligns cds_pairs 2>/dev/null || echo "No intermediates to tar or tar failed"
# optionally remove intermediate dirs to save node local space
rm -rf protein_aligns codon_aligns cds_pairs || true

# ========== 7) copy important files back to RESULTS_ROOT (done in trap as well) ==========
echo "[7] copying important outputs to $RESULTS_ROOT/$id/"
mkdir -p "${RESULTS_ROOT}/${id}"

# copy BLAST outputs
cp -a hap1_vs_hap2.blast hap2_vs_hap1.blast hap1.top.tsv hap2.top.tsv "${RESULTS_ROOT}/${id}/" 2>/dev/null || true

# copy PAV outputs
if [[ -d pav_blast_out ]]; then
    rsync -a pav_blast_out "${RESULTS_ROOT}/${id}/" || cp -a pav_blast_out "${RESULTS_ROOT}/${id}/"
fi

# copy merged axt and kaks
if [[ -s "${id}.merged.axt" ]]; then
    cp -a "${id}.merged.axt" "${RESULTS_ROOT}/${id}/" 2>/dev/null || true
fi
if [[ -s "${id}.kaks.csv" ]]; then
    cp -a "${id}.kaks.csv" "${RESULTS_ROOT}/${id}/" 2>/dev/null || true
fi

# copy intermediate tarball
if [[ -s "${INTERMEDIATE_TAR}" ]]; then
    cp -a "${INTERMEDIATE_TAR}" "${RESULTS_ROOT}/${id}/" 2>/dev/null || true
fi
# remove WORKDIR to free node local space (comment out if you want to inspect)
    rm -rf "${WORKDIR}" || true
    echo "[$(date)] trap complete."
echo "[$(date)] DONE sample $id"
# trap will also try to copy back/cleanup
exit 0
