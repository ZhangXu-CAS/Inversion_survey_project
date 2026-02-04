#!/bin/bash
#SBATCH --job-name=1031braker1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=60G
#SBATCH --time=12:00:00
#SBATCH --account=def-rieseber
#SBATCH --array=0-48
#SBATCH --out=log/braker1_min%j_%a.log

set -euo pipefail
# set -x   # uncomment for debugging single job

module load apptainer/1.3.5

# ---- Config ----
THREADS="${SLURM_CPUS_PER_TASK:-1}"
PROJ_DIR="/project/6004758/zhangxu/TE_anno"
cd "$PROJ_DIR"

# sample id
id=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" braker1_ids)
echo "Start sample: $id  JOB=$SLURM_JOB_ID TASK=$SLURM_ARRAY_TASK_ID HOST=$(hostname) TIME=$(date)"

if [[ -f "${id}/${id}_hap1_revchr.fa.masked" ]]; then
    GENOME_MASK="${id}/${id}_hap1_revchr.fa.masked"
else
    GENOME_MASK="../4Annotation_done/${id}/${id}_hap1_revchr.fa"
fi

OUT_SUBDIR="${id}/braker_out1"
mkdir -p "$OUT_SUBDIR"

# node-local tmp dir on host (created here)
TMP_BASE="${SLURM_TMPDIR:-/tmp}"
TMP_WORK=$(mktemp -d "${TMP_BASE}/braker_${id}_XXXXXXXX")
echo "TMP_WORK (host) = $TMP_WORK"

# container and binds
CONTAINER_IMG="/project/6004758/zhangxu/4Annotation_done/braker3.sif"
AUG_CFG_HOST="${HOME}/augustus_config_copy/config"
GENEMARK_BIN_HOST="/home/zhangxu/software/GeneMark-ETP/bin"

# basic input checks
if [[ ! -f "$GENOME_MASK" ]]; then
  echo "ERROR: genome not found: $GENOME_MASK" >&2
  exit 2
fi
if [[ ! -f "$CONTAINER_IMG" ]]; then
  echo "ERROR: container image not found: $CONTAINER_IMG" >&2
  exit 2
fi
if [[ ! -d "$AUG_CFG_HOST" ]]; then
  echo "ERROR: AUGUSTUS config host dir not found: $AUG_CFG_HOST" >&2
  exit 2
fi
if [[ ! -d "$GENEMARK_BIN_HOST" ]]; then
  echo "ERROR: GeneMark bin host dir not found: $GENEMARK_BIN_HOST" >&2
  exit 2
fi

KEY_FILES=(braker.gtf braker.codingseq braker.aa braker.log augustus.hints.gff)

copy_back_minimal() {
  rc=$?
  echo "[$(date)] exit handler (rc=$rc): copying minimal key files back..."
  mkdir -p "$OUT_SUBDIR"

  STAGE_DIR=$(mktemp -d "${TMP_BASE}/braker_stage_${id}_XXXXXXXX")
  echo "Staging in $STAGE_DIR"

  for f in "${KEY_FILES[@]}"; do
    if [[ -f "${TMP_WORK}/${f}" ]]; then
      cp -v "${TMP_WORK}/${f}" "${STAGE_DIR}/"
    else
      echo "Note: ${f} not found in TMP_WORK"
    fi
  done

  TAR_NAME="${OUT_SUBDIR}/${id}_braker_minimal_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.tar.gz"
  if [[ -n "$(ls -A "$STAGE_DIR" 2>/dev/null || true)" ]]; then
    tar -C "$STAGE_DIR" -czf "$TAR_NAME" . || echo "Warning: tar failed"
    echo "Created minimal tar: $TAR_NAME  size: $(du -h "$TAR_NAME" 2>/dev/null | cut -f1 || echo 'N/A')"
    for f in braker.gtf braker.codingseq braker.aa braker.log; do
      if [[ -f "${STAGE_DIR}/${f}" ]]; then
        cp -v "${STAGE_DIR}/${f}" "${OUT_SUBDIR}/"
      fi
    done
  else
    echo "No key files found to copy back."
  fi

  rm -rf "$STAGE_DIR"
  if [[ -d "$TMP_WORK" ]]; then
    rm -rf "$TMP_WORK" || echo "Warning: failed to remove $TMP_WORK"
  fi

  echo "Exit handler finished."
  return $rc
}
trap copy_back_minimal EXIT

echo "Running BRAKER in node-local TMP ($TMP_WORK) with threads=$THREADS..."

# bind the host TMP_WORK into the container at a safe container-internal path
CONTAINER_TMP="/tmp/braker_${id}"

apptainer exec \
  --bind "${PROJ_DIR}:/work:rw" \
  --bind "${AUG_CFG_HOST}:/augustus_config:ro" \
  --bind "${GENEMARK_BIN_HOST}:/opt/genemark:ro" \
  --bind "${TMP_WORK}:${CONTAINER_TMP}:rw" \
  "$CONTAINER_IMG" \
  bash -c "
    set -euo pipefail
    cd /work
    # NOTE: the bound CONTAINER_TMP already exists and is writable
    /usr/bin/env braker.pl --genome=\"/work/${GENOME_MASK}\" \
      --species=\"${id}_hap1\" \
      --threads=${THREADS} \
      --useexisting \
      --GENEMARK_PATH=/opt/genemark/gmes \
      --AUGUSTUS_CONFIG_PATH=/augustus_config \
      --workingdir=\"${CONTAINER_TMP}\"
  "

echo "BRAKER run finished for $id at $(date)"
exit 0
