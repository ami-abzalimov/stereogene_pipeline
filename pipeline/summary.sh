#!/usr/bin/env bash
set -euo pipefail
trap 'log_err "[ERROR] Failure at line $LINENO"; exit 1' ERR

export OMP_NUM_THREADS=1

# Config vars
CSV="encode_filess.csv"
API="https://www.encodeproject.org"
RAW_CHR_SIZES="hg38.chrom.sizes"
CHR_SIZES="$RAW_CHR_SIZES"
OUTROOT="encode_output"
MAX_RETRIES=3
RETRY_DELAY=5

# Binary paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_BINNER="$SCRIPT_DIR/src/Binner"
BIN_SMOOTHER="$SCRIPT_DIR/src/Smoother"
BIN_CONFOUND="$SCRIPT_DIR/src/Confounder"
BIN_PROJECT="$SCRIPT_DIR/src/Projector"
BIN_SG="$SCRIPT_DIR/src/StereoGene"
SG_OPTS="outRes=TAB outDistr=1 outChrom=1 outPrjBGr=1 Distances=1 plotType=NONE CP=0"

# Logging helpers
log_ts(){ date '+%Y-%m-%d %H:%M:%S'; }
log_info(){ printf '%s [INFO]  %s\n' "$(log_ts)" "$*"; }
log_warn(){ printf '%s [WARN]  %s\n' "$(log_ts)" "$*"; }
log_err(){ printf '%s [ERROR] %s\n' "$(log_ts)" "$*"; }

# Minimal preflight
if [[ ! -f "$CSV" ]]; then
  log_err "CSV file not found: $CSV"
  exit 1
fi
if ! command -v aria2c >/dev/null 2>&1; then
  log_err "aria2c required but not found"
  exit 1
fi

# Chrom sizes
if [[ ! -f "$CHR_SIZES" ]]; then
  log_warn "Chrom sizes missing: $CHR_SIZES"
fi

mkdir -p "$OUTROOT"
PROCESSED_FILE="$OUTROOT/.processed_ids"
: > "$PROCESSED_FILE"

# Filename sanitize
safe_basename(){
  local name="$1"
  name="${name##*/}"
  name="${name%%\?*}"
  name="${name//[$'\r\n\t ']/_}"
  name="$(printf '%s' "$name" | sed 's/[^A-Za-z0-9._-]/_/g')"
  [[ -z "$name" ]] && name="file_$(date +%s%N)"
  printf '%s' "$name"
}

# Mark processed
mark_processed(){
  local id="$1"
  mkdir -p "$(dirname "$PROCESSED_FILE")"
  grep -qxF "$id" "$PROCESSED_FILE" 2>/dev/null || printf '%s\n' "$id" >> "$PROCESSED_FILE"
}
# Check processed
already_processed(){
  local id="$1"
  [[ -f "$PROCESSED_FILE" ]] && grep -qxF "$id" "$PROCESSED_FILE" 2>/dev/null
}

# Simple download
download_with_aria2(){
  local id="$1" dir="$2" url="$API/files/$id/@@download"
  mkdir -p "$dir"
  log_info "Downloading $id -> $dir"
  aria2c --file-allocation=none -x16 -s16 -k1M --max-tries="$MAX_RETRIES" --retry-wait="$RETRY_DELAY" --dir="$dir" "$url" >/dev/null 2>&1 || true
  find "$dir" -maxdepth 1 -type f -printf '%T@ %p\n' 2>/dev/null | sort -nr | head -n1 | cut -d' ' -f2- || true
}

# Track processing
process_track(){
  local raw="$1" tgt_dir="$2"
  raw="${raw#"${raw%%[![:space:]]*}"}"
  raw="${raw%"${raw##*[![:space:]]}"}"
  raw="${raw#*://}"
  raw="${raw%%\?*}"
  raw="${raw%/}"
  raw="${raw%,}"
  local id
  id=$(printf '%s' "$raw" | grep -oE 'ENCFF[0-9A-Z]+' || true)
  if [[ -z "$id" ]]; then id="$(basename "$raw")"; fi
  id=$(echo "$id" | sed 's/\?.*$//' | sed 's#/.*$##' | sed 's/\.[^./?]*$//' | tr '[:lower:]' '[:upper:]')
  log_info "Process ID $id -> $tgt_dir"
  if already_processed "$id"; then log_info "$id skipped"; return 0; fi

  mkdir -p "$tgt_dir"
  local dl
  dl=$(download_with_aria2 "$id" "$tgt_dir" || true)
  if [[ -z "$dl" || ! -f "$dl" ]]; then
    log_warn "No download for $id"
    mark_processed "$id"
    return 0
  fi
  local dlfile="$dl"
  log_info "Downloaded file: $dlfile"

  local base fmt name_noext inner final_wig=""
  base=$(basename "$dlfile")
  fmt="${base##*.}"
  name_noext="${base%.*}"

  case "${fmt,,}" in
    bam)
      # Convert BAM
      log_info "BAM -> wig"
      samtools sort -@4 -o "${tgt_dir}/${name_noext}_sorted.bam" "$dlfile"
      samtools index "${tgt_dir}/${name_noext}_sorted.bam"
      bamToBed -i "${tgt_dir}/${name_noext}_sorted.bam" > "${tgt_dir}/${name_noext}.bed"
      genomeCoverageBed -bg -i "${tgt_dir}/${name_noext}.bed" -g "$CHR_SIZES" > "${tgt_dir}/${name_noext}.bedGraph"
      sort -k1,1 -k2,2n "${tgt_dir}/${name_noext}.bedGraph" > "${tgt_dir}/${name_noext}_sorted.bedGraph"
      bedGraphToBigWig "${tgt_dir}/${name_noext}_sorted.bedGraph" "$CHR_SIZES" "${tgt_dir}/${name_noext}.bw"
      bigWigToWig "${tgt_dir}/${name_noext}.bw" "${tgt_dir}/${name_noext}.wig"
      rm -f "${tgt_dir}/${name_noext}.bed" "${tgt_dir}/${name_noext}.bedGraph" "${tgt_dir}/${name_noext}_sorted.bedGraph" "${tgt_dir}/${name_noext}.bw" "${tgt_dir}/${name_noext}_sorted.bam" "${tgt_dir}/${name_noext}_sorted.bam.bai"
      final_wig="${tgt_dir}/${name_noext}.wig"
      ;;
    bed)
      # Convert BED
      log_info "BED -> wig"
      genomeCoverageBed -bg -i "$dlfile" -g "$CHR_SIZES" > "${tgt_dir}/${name_noext}.bedGraph"
      sort -k1,1 -k2,2n "${tgt_dir}/${name_noext}.bedGraph" > "${tgt_dir}/${name_noext}_sorted.bedGraph"
      bedGraphToBigWig "${tgt_dir}/${name_noext}_sorted.bedGraph" "$CHR_SIZES" "${tgt_dir}/${name_noext}.bw"
      bigWigToWig "${tgt_dir}/${name_noext}.bw" "${tgt_dir}/${name_noext}.wig"
      rm -f "${tgt_dir}/${name_noext}.bedGraph" "${tgt_dir}/${name_noext}_sorted.bedGraph" "${tgt_dir}/${name_noext}.bw"
      final_wig="${tgt_dir}/${name_noext}.wig"
      ;;
    bb|bigbed)
      # Convert BigBed
      log_info "BigBed -> wig"
      bigBedToBed "$dlfile" "${tgt_dir}/${name_noext}.bed"
      genomeCoverageBed -bg -i "${tgt_dir}/${name_noext}.bed" -g "$CHR_SIZES" > "${tgt_dir}/${name_noext}.bedGraph"
      sort -k1,1 -k2,2n "${tgt_dir}/${name_noext}.bedGraph" > "${tgt_dir}/${name_noext}_sorted.bedGraph"
      bedGraphToBigWig "${tgt_dir}/${name_noext}_sorted.bedGraph" "$CHR_SIZES" "${tgt_dir}/${name_noext}.bw"
      bigWigToWig "${tgt_dir}/${name_noext}.bw" "${tgt_dir}/${name_noext}.wig"
      rm -f "${tgt_dir}/${name_noext}.bed" "${tgt_dir}/${name_noext}.bedGraph" "${tgt_dir}/${name_noext}_sorted.bedGraph" "${tgt_dir}/${name_noext}.bw"
      final_wig="${tgt_dir}/${name_noext}.wig"
      ;;
    bw|bigwig)
      # Convert BigWig
      log_info "BigWig -> wig"
      bigWigToWig "$dlfile" "${tgt_dir}/${name_noext}.wig"
      final_wig="${tgt_dir}/${name_noext}.wig"
      ;;
    gz)
      # Handle gzip
      inner="${name_noext}"
      case "${inner##*.}" in
        bed)
          log_info "bed.gz -> wig"
          gunzip -c "$dlfile" > "${tgt_dir}/${inner}"
          genomeCoverageBed -bg -i "${tgt_dir}/${inner}" -g "$CHR_SIZES" > "${tgt_dir}/${inner}.bedGraph"
          sort -k1,1 -k2,2n "${tgt_dir}/${inner}.bedGraph" > "${tgt_dir}/${inner}_sorted.bedGraph"
          bedGraphToBigWig "${tgt_dir}/${inner}_sorted.bedGraph" "$CHR_SIZES" "${tgt_dir}/${inner}.bw"
          bigWigToWig "${tgt_dir}/${inner}.bw" "${tgt_dir}/${inner}.wig"
          rm -f "${tgt_dir}/${inner}.bed" "${tgt_dir}/${inner}.bedGraph" "${tgt_dir}/${inner}_sorted.bedGraph" "${tgt_dir}/${inner}.bw"
          final_wig="${tgt_dir}/${inner}.wig"
          ;;
        wig)
          log_info "wig.gz -> wig"
          gunzip -c "$dlfile" > "${tgt_dir}/${inner}"
          final_wig="${tgt_dir}/${inner}"
          ;;
        *)
          log_warn "Unknown gz inner: $dlfile"
          final_wig=""
          ;;
      esac
      ;;
    *)
      # Unknown extension
      log_warn "Unhandled ext: $dlfile"
      final_wig=""
      ;;
  esac

  if [[ -n "${final_wig}" && -f "$final_wig" ]]; then
    # Postprocess wig
    [[ -x "$BIN_BINNER" ]] && "$BIN_BINNER" "$final_wig" || log_warn "Binner absent"
    [[ -x "$BIN_SMOOTHER" ]] && "$BIN_SMOOTHER" "$final_wig" || log_warn "Smoother absent"
    mark_processed "$id"
    printf '%s' "$final_wig"
    return 0
  else
    log_warn "No wig produced for $id"
    mark_processed "$id"
    return 0
  fi
}

# StereoGene wrapper
run_sg(){
  local f1="$1" f2="$2" outdir="$3"
  mkdir -p "$outdir"
  log_info "StereoGene: $(basename "$f1") vs $(basename "$f2")"
  [[ -x "$BIN_SG" ]] && "$BIN_SG" chrom="$CHR_SIZES" $SG_OPTS "$f1" "$f2" || log_warn "StereoGene absent"
}

# Main CSV processing
log_info "Start CSV processing"
while IFS=$'\t' read -r assay target files || [[ -n "$assay" ]]; do
  [[ -z "$assay" ]] && continue
  [[ -z "$target" ]] && continue
  log_info "Row: $assay | $target"
  base_dir="$OUTROOT/${target// /_}"
  assay_dir="$base_dir/${assay// /_}"
  mkdir -p "$assay_dir"
  IFS=',' read -r -a arr <<< "$files"
  row_wigs=()
  for raw in "${arr[@]}"; do
    raw="${raw#"${raw%%[![:space:]]*}"}"
    raw="${raw%"${raw##*[![:space:]]}"}"
    [[ -z "$raw" ]] && continue
    wigpath=$(process_track "$raw" "$assay_dir" || true)
    if [[ -n "$wigpath" && -f "$wigpath" ]]; then
      log_info "Got wig: $wigpath"
      row_wigs+=( "$wigpath" )
    fi
  done

  if [[ ${#row_wigs[@]} -gt 0 ]]; then
    tmp_list="$(mktemp -t list.wig.XXXXXX)"
    printf '%s\n' "${row_wigs[@]}" > "$tmp_list"
    log_info "Run Confounder/Projector"
    [[ -x "$BIN_CONFOUND" ]] && "$BIN_CONFOUND" "$tmp_list" || log_warn "Confounder absent"
    [[ -x "$BIN_PROJECT" ]] && "$BIN_PROJECT" "$assay_dir" || log_warn "Projector absent"
    rm -f "$tmp_list"
  else
    log_info "No wigs -> skip confound/project"
  fi
done < <(tail -n +2 "$CSV")

# StereoGene pairwise
log_info "StereoGene pairwise"
for base_dir in "$OUTROOT"/*; do
  [[ -d "$base_dir" ]] || continue
  for assay_dir in "$base_dir"/*; do
    [[ -d "$assay_dir" ]] || continue
    mapfile -t wigf < <(find "$assay_dir" -maxdepth 1 -name '*.wig' -print 2>/dev/null || true)
    [[ ${#wigf[@]} -lt 2 ]] && continue
    log_info "Pairwise in $assay_dir"
    for ((i=0;i<${#wigf[@]};i++)); do
      for ((j=i+1;j<${#wigf[@]};j++)); do
        run_sg "${wigf[i]}" "${wigf[j]}" "$assay_dir/sg_results" || log_warn "StereoGene pair failed"
      done
    done
  done
done

# Archive assays
log_info "Archive assays"
for base_dir in "$OUTROOT"/*; do
  [[ -d "$base_dir" ]] || continue
  for assay_dir in "$base_dir"/*; do
    [[ -d "$assay_dir" ]] || continue
    if find "$assay_dir" -mindepth 1 -print -quit >/dev/null 2>&1; then
      tarball="$assay_dir.tar.gz"
      log_info "Create $tarball"
      tar -czf "$tarball" -C "$base_dir" "$(basename "$assay_dir")" || log_warn "Tar failed"
    else
      log_info "Empty dir skip"
    fi
  done
done

log_info "Pipeline done"
