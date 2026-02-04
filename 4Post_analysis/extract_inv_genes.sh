#!/usr/bin/env bash
set -euo pipefail
module load bedtools
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 All_PASS_inversions.bed braker_dir outdir"
  exit 1
fi

INV_BED="$1"
BRAKER_DIR="$2"
OUTDIR="$3"
mkdir -p "${OUTDIR}"

TMPDIR=$(mktemp -d)
cleanup(){ rm -rf "${TMPDIR}"; }
trap cleanup EXIT

debug_log="${OUTDIR}/extract_inv_genes.debug.log"
: > "${debug_log}"

echo "Start: $(date)" | tee -a "${debug_log}"
echo "Inputs: inv_bed=${INV_BED}, braker_dir=${BRAKER_DIR}, outdir=${OUTDIR}" | tee -a "${debug_log}"

# --- auto-detect 0-based vs 1-based for inversion file ---
# get a sample of start coords and check if any zero exists
min_start=$(awk 'NR>0{if($2+0<min||NR==1)min=$2+0}END{print min}' "${INV_BED}" 2>/dev/null || echo "")
if [ -z "${min_start}" ]; then
  echo "ERROR: could not read ${INV_BED} or it's empty." | tee -a "${debug_log}"
  exit 1
fi

if [ "${min_start}" -eq 0 ]; then
  INV_IS_0BASED=true
  echo "Detected inversion BED appears 0-based (min start == 0)." | tee -a "${debug_log}"
else
  INV_IS_0BASED=false
  echo "Detected inversion BED appears 1-based (min start > 0). Will convert to 0-based for bedtools." | tee -a "${debug_log}"
fi

# prepare sample list (col4)
cut -f4 "${INV_BED}" | sort -u > "${TMPDIR}/sample_ids.txt"

# function: extract gene spans robustly; writes out_bed (chr\tstart0\tend\tgeneid)
extract_gene_spans() {
  local gtf="$1"
  local out_bed="$2"
  local tmp="${out_bed}.tmp"

  : > "${tmp}"

  # 1) try feature gene with common attribute styles
  awk 'BEGIN{FS="\t"; OFS="\t"}
    $0 ~ /^#/ { next }
    $3=="gene" {
      attr=$9
      gid=""
      if (match(attr, /gene_id[ \t]*"([^"]+)"/, a)) gid=a[1]
      else if(match(attr, /gene_id=([^; \t]+)/, a)) gid=a[1]
      else if(match(attr, /ID=([^; \t]+)/, a)) gid=a[1]
      if(gid!=""){
        start=$4-1; end=$5+0
        printf("%s\t%d\t%d\t%s\n", $1, start, end, gid) >> "'"${tmp}"'"
      }
    }' "${gtf}"

  if [ -s "${tmp}" ]; then
    mv "${tmp}" "${out_bed}"
    return 0
  fi
  rm -f "${tmp}"

  # 2) try mRNA/transcript and aggregate to gene (by gene_id or ID)
  awk 'BEGIN{FS="\t"; OFS="\t"}
    $0 ~ /^#/ { next }
    ($3=="mRNA" || $3=="transcript" || $3=="tRNA" || $3=="mRNA_template") {
      attr=$9; gid=""
      if (match(attr, /gene_id[ \t]*"([^"]+)"/, a)) gid=a[1]
      else if(match(attr, /gene_id=([^; \t]+)/, a)) gid=a[1]
      else if(match(attr, /Parent=([^; \t]+)/, a)) gid=a[1]
      else if(match(attr, /ID=([^; \t]+)/, a)) gid=a[1]
      if(gid!=""){ start=$4-1; end=$5+0; print $1, start, end, gid >> "'"${tmp}"'" }
    }' "${gtf}"

  if [ -s "${tmp}" ]; then
    # aggregate per gene id -> min start, max end
    sort -k1,1 -k4,4 -k2,2n "${tmp}" \
      | awk 'BEGIN{FS=OFS="\t"}{
          key=$1 FS $4
          if(key!=prev){
            if(NR>1) print ochr, ost, oend, og;
            ochr=$1; ost=$2; oend=$3; og=$4; prev=key;
          } else {
            if($2 < ost) ost=$2;
            if($3 > oend) oend=$3;
          }
        } END{ if(prev!="") print ochr, ost, oend, og }' > "${out_bed}"
    rm -f "${tmp}"
    return 0
  fi
  rm -f "${tmp}"

  # 3) fallback CDS
  awk 'BEGIN{FS="\t"; OFS="\t"}
    $0 ~ /^#/ { next }
    $3=="CDS" {
      attr=$9; gid=""
      if (match(attr, /gene_id[ \t]*"([^"]+)"/, a)) gid=a[1]
      else if(match(attr, /gene_id=([^; \t]+)/, a)) gid=a[1]
      else if(match(attr, /Parent=([^; \t]+)/, a)) gid=a[1]
      else if(match(attr, /ID=([^; \t]+)/, a)) gid=a[1]
      if(gid!=""){ start=$4-1; end=$5+0; print $1, start, end, gid >> "'"${tmp}"'" }
    }' "${gtf}"

  if [ -s "${tmp}" ]; then
    sort -k1,1 -k4,4 -k2,2n "${tmp}" \
      | awk 'BEGIN{FS=OFS="\t"}{
          key=$1 FS $4
          if(key!=prev){
            if(NR>1) print ochr, ost, oend, og;
            ochr=$1; ost=$2; oend=$3; og=$4; prev=key;
          } else {
            if($2 < ost) ost=$2;
            if($3 > oend) oend=$3;
          }
        } END{ if(prev!="") print ochr, ost, oend, og }' > "${out_bed}"
    rm -f "${tmp}"
    return 0
  fi
  rm -f "${tmp}"
  # nothing found
  touch "${out_bed}"
  return 0
}

# helper check: validate bed-like file has integer cols 2 and 3 and is tab-separated
validate_bed() {
  local f="$1"
  if [ ! -s "${f}" ]; then
    return 1
  fi
  # check tab present on first line and numeric cols
  head -n 5 "${f}" | awk -F"\t" '{
    if(NF < 4){ exit 2 }
    if($2+0 != $2 || $3+0 != $3){ exit 3 }
  }' || return $?
  return 0
}

# iterate samples
while read -r ID; do
  [ -z "$ID" ] && continue
  echo "Processing sample: ${ID}" | tee -a "${debug_log}"

  samp_inv="${TMPDIR}/${ID}.inv.raw.bed"
  awk -v id="$ID" '$4==id{print $1"\t"$2"\t"$3"\t"$4}' "${INV_BED}" > "${samp_inv}" || true

  if [ ! -s "${samp_inv}" ]; then
    echo "  WARN: no inversions for ${ID}" | tee -a "${debug_log}"
    continue
  fi

  # convert inv to 0-based bed for bedtools if needed (auto-detect)
  samp_inv_bed="${TMPDIR}/${ID}.inv.bed"
  if [ "${INV_IS_0BASED}" = true ]; then
    cp "${samp_inv}" "${samp_inv_bed}"
  else
    awk 'BEGIN{OFS="\t"} { s=$2-1; if(s<0) s=0; print $1, s, $3, $4 }' "${samp_inv}" > "${samp_inv_bed}"
  fi

  # GTF path
  gtf="${BRAKER_DIR}/${ID}.braker.gtf"
  if [ ! -s "${gtf}" ]; then
    echo "  WARN: missing GTF: ${gtf}" | tee -a "${debug_log}"
    continue
  fi

  genes_bed="${TMPDIR}/${ID}.genes.bed"
  extract_gene_spans "${gtf}" "${genes_bed}"

  # validate genes_bed
  if ! validate_bed "${genes_bed}"; then
    echo "  ERROR: generated genes_bed invalid or empty for ${ID}. See debug log." | tee -a "${debug_log}"
    echo "  Path: ${genes_bed}" | tee -a "${debug_log}"
    echo "  First 20 lines of GTF (for quick check):" >> "${debug_log}"
    sed -n '1,20p' "${gtf}" >> "${debug_log}"
    echo "  Dump of generated genes_bed (first 20 lines):" >> "${debug_log}"
    sed -n '1,20p' "${genes_bed}" >> "${debug_log}"
    echo "  Suggestion: inspect attribute formats (gene_id, ID=, gene_id=) in GTF. Example grep:" >> "${debug_log}"
    echo "    grep -m 5 -E 'gene_id|ID=' ${gtf} | sed -n '1,10p'" >> "${debug_log}"
    continue
  fi

  # find overlapping gene IDs with inversions
  overlap_ids="${TMPDIR}/${ID}.inv_gene_ids.txt"
  bedtools intersect -a "${genes_bed}" -b "${samp_inv_bed}" -wa -u 2>>"${debug_log}" \
    | cut -f4 | sort -u > "${overlap_ids}" || true

  out_gtf="${OUTDIR}/${ID}.invgene.gtf"
  : > "${out_gtf}"

  if [ -s "${overlap_ids}" ]; then
    # extract full GTF lines that match gene ids in overlap_ids
    awk 'BEGIN{ while(getline<ARGV[2]) ids[$1]=1; close(ARGV[2]) }
         {
           if(match($0, /gene_id[ \t]*"([^"]+)"/, a) && a[1] in ids) print $0
           else if(match($0, /gene_id=([^; \t]+)/, b) && b[1] in ids) print $0
           else if(match($0, /ID=([^; \t]+)/, c) && c[1] in ids) print $0
         }' "${gtf}" "${overlap_ids}" > "${out_gtf}"
    echo "  wrote ${out_gtf} (genes overlapping inversions)" | tee -a "${debug_log}"
  else
    echo "  No genes overlap inversions for ${ID}." | tee -a "${debug_log}"
    : > "${out_gtf}"
  fi

  # downstream right 1kb of inversion start (convert start back to original start location if we subtracted)
  right_flank_bed="${TMPDIR}/${ID}.inv_left1kb.bed"
  awk 'BEGIN{OFS="\t"} { s=$2+1000; if(s<0) s=0; e=$2-1; if(e>=0) print $1, s, e, $4 }' "${samp_inv_bed}" > "${right_flank_bed}"

  up_ids="${TMPDIR}/${ID}.inv_up_gene_ids.txt"
  bedtools intersect -a "${genes_bed}" -b "${right_flank_bed}" -wa -u 2>>"${debug_log}" \
    | cut -f4 | sort -u > "${up_ids}" || true

  out_down_gtf="${OUTDIR}/${ID}.inv_downgene.gtf"
  : > "${out_down_gtf}"
  if [ -s "${up_ids}" ]; then
    awk 'BEGIN{ while(getline<ARGV[2]) ids[$1]=1; close(ARGV[2]) }
         {
           if(match($0, /gene_id[ \t]*"([^"]+)"/, a) && a[1] in ids) print $0
           else if(match($0, /gene_id=([^; \t]+)/, b) && b[1] in ids) print $0
           else if(match($0, /ID=([^; \t]+)/, c) && c[1] in ids) print $0
         }' "${gtf}" "${up_ids}" > "${out_down_gtf}"
    echo "  wrote ${out_down_gtf} (genes in right 1kb downstream)" | tee -a "${debug_log}"
  else
    echo "  No genes found in 1kb downstream region for ${ID}." | tee -a "${debug_log}"
    : > "${out_down_gtf}"
  fi

done < "${TMPDIR}/sample_ids.txt"

echo "Done: $(date)" | tee -a "${debug_log}"
echo "If any samples produced errors, check ${debug_log} for details."
