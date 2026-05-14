#!/usr/bin/env bash
# organize_rds.sh
#
# Migrates scDiffCom .rds files from source to per-dataset destination
# directories, then audits for genes missing from any of the three datasets.
#
# Usage:
#   bash organize_rds.sh                   # live run
#   DRY_RUN=true bash organize_rds.sh      # preview only — no files moved
#   REPORT_FILE=missing.txt bash organize_rds.sh   # also write audit to file

set -euo pipefail

# ─── Configuration ────────────────────────────────────────────────────────────
DRY_RUN="${DRY_RUN:-false}"
REPORT_FILE="${REPORT_FILE:-}"

SOURCE_DIR="${HOME}/CCC-scDiffCom/results"
TARGET_ROOT="${HOME}/Thesis/CCC/outputs/scDiffCom/scDiffComs/Consensus_Cell_Type"

DATASETS=(
    "Bassez2021_Breast_step4"
    "Qian2020_Breast_step4"
    "Wu2021_Breast_step4"
)

SUFFIX="_scDiffCom.rds"

# ─── Helpers ──────────────────────────────────────────────────────────────────
log()  { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*"; }
info() { log "INFO  $*"; }
warn() { log "WARN  $*" >&2; }

# Print command in dry-run mode; execute it otherwise.
run() {
    if [[ "$DRY_RUN" == "true" ]]; then
        log "DRY   $(printf '%q ' "$@")"
    else
        "$@"
    fi
}

# ─── Pre-flight checks ────────────────────────────────────────────────────────
if [[ ! -d "$SOURCE_DIR" ]]; then
    warn "Source directory not found: $SOURCE_DIR"
    exit 1
fi

if ! command -v rsync &>/dev/null; then
    warn "rsync not found — install it or replace with cp"
    exit 1
fi

[[ "$DRY_RUN" == "true" ]] && info "*** DRY RUN MODE — no files will be moved ***"

# ─── Phase 1: File Migration ───────────────────────────────────────────────────
info "=== Phase 1: File Migration ==="
info "  Source : $SOURCE_DIR"
info "  Target : $TARGET_ROOT"
echo

moved=0; skipped=0; unmatched=0; errors=0

while IFS= read -r -d '' src_file; do
    filename="$(basename "$src_file")"

    # Match filename to a known dataset by checking the known suffix patterns.
    # Gene names may contain underscores, so we anchor on the *known* dataset
    # strings rather than splitting naively on '_'.
    matched_dataset=""
    for ds in "${DATASETS[@]}"; do
        if [[ "$filename" == *"_${ds}${SUFFIX}" ]]; then
            matched_dataset="$ds"
            break
        fi
    done

    if [[ -z "$matched_dataset" ]]; then
        warn "No dataset match for '$filename' — skipping"
        (( unmatched++ )) || true
        continue
    fi

    dest_dir="${TARGET_ROOT}/${matched_dataset}"
    dest_file="${dest_dir}/${filename}"

    # Create destination directory when it doesn't exist yet.
    if [[ "$DRY_RUN" == "true" ]]; then
        [[ -d "$dest_dir" ]] || log "DRY   mkdir -p '$dest_dir'"
    else
        mkdir -p "$dest_dir"
    fi

    # Timestamp-based override: skip if destination exists and is not older.
    # -nt returns true when src is newer than dest (mtime comparison).
    if [[ -e "$dest_file" && ! "$src_file" -nt "$dest_file" ]]; then
        info "SKIP   (dest is current or newer) : $filename"
        (( skipped++ )) || true
        continue
    fi

    if [[ -e "$dest_file" ]]; then
        info "UPDATE (src is newer)             : $filename  →  ${matched_dataset}/"
    else
        info "COPY   (new file)                 : $filename  →  ${matched_dataset}/"
    fi

    # rsync -a preserves timestamps/permissions; --checksum skips byte-identical
    # files even if timestamps differ (belt-and-suspenders check).
    if [[ "$DRY_RUN" == "true" ]]; then
        log "DRY   rsync --archive --checksum '$src_file' '$dest_file'"
        (( moved++ )) || true
    else
        if rsync --archive --checksum "$src_file" "$dest_file"; then
            (( moved++ )) || true
        else
            warn "rsync failed for: $filename"
            (( errors++ )) || true
        fi
    fi

done < <(find "$SOURCE_DIR" -maxdepth 1 -name "*${SUFFIX}" -print0 | sort -z)

echo
info "Migration summary"
info "  Moved/updated : ${moved}"
info "  Skipped       : ${skipped}  (dest already current)"
info "  Unmatched     : ${unmatched}  (no dataset found in filename)"
info "  Errors        : ${errors}"
echo

# ─── Phase 2: Missingness Audit ───────────────────────────────────────────────
# Scans the *destination* directories so the audit reflects current on-disk
# state (including files placed there by previous runs).

info "=== Phase 2: Missingness Audit (scanning destination) ==="
echo

declare -A gene_present   # gene → newline-separated datasets already seen

for ds in "${DATASETS[@]}"; do
    ds_dir="${TARGET_ROOT}/${ds}"

    if [[ ! -d "$ds_dir" ]]; then
        warn "  Dataset dir not found — skipping audit for: $ds_dir"
        continue
    fi

    file_count=0
    while IFS= read -r -d '' f; do
        fname="$(basename "$f")"
        # Strip the known trailing pattern to isolate the gene name.
        gene="${fname%_${ds}${SUFFIX}}"

        if [[ -z "$gene" ]]; then
            warn "  Could not extract gene name from: $fname"
            continue
        fi

        (( file_count++ )) || true

        if [[ -v gene_present["$gene"] ]]; then
            gene_present["$gene"]+=$'\n'"$ds"
        else
            gene_present["$gene"]="$ds"
        fi
    done < <(find "$ds_dir" -maxdepth 1 -name "*_${ds}${SUFFIX}" -print0)

    info "  ${ds}: ${file_count} file(s) found"
done

echo

# ─── Build report ─────────────────────────────────────────────────────────────
SEP="━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
mapfile -t report_lines < <(
    echo "$SEP"
    printf "  MISSINGNESS REPORT  (generated %s)\n" "$(date '+%Y-%m-%d %H:%M:%S')"
    echo "$SEP"
    printf "  %-40s  %s\n" "GENE" "MISSING FROM"
    echo "  ──────────────────────────────────────────────────────────────────────"
)

total_genes=0
complete_genes=0
incomplete_genes=0
incomplete_lines=()

for gene in $(printf '%s\n' "${!gene_present[@]}" | sort); do
    (( total_genes++ )) || true

    present_str="${gene_present[$gene]}"

    missing_list=()
    for ds in "${DATASETS[@]}"; do
        # Check whether this dataset appears as a whole entry in the record.
        found=false
        while IFS= read -r entry; do
            [[ "$entry" == "$ds" ]] && found=true && break
        done <<< "$present_str"
        $found || missing_list+=("$ds")
    done

    if (( ${#missing_list[@]} == 0 )); then
        (( complete_genes++ )) || true
    else
        (( incomplete_genes++ )) || true
        missing_str="$(IFS=', '; echo "${missing_list[*]}")"
        incomplete_lines+=("$(printf '  %-40s  %s' "$gene" "$missing_str")")
    fi
done

# Append incomplete-gene lines (or a "none" message) then the summary.
if (( ${#incomplete_lines[@]} == 0 )); then
    report_lines+=("  (none — all genes are complete across all datasets)")
else
    report_lines+=("${incomplete_lines[@]}")
fi

report_lines+=(
    "$SEP"
    "$(printf '  %-20s %d' 'Total genes:'   "$total_genes")"
    "$(printf '  %-20s %d  (present in all %d datasets)' 'Complete:' "$complete_genes" "${#DATASETS[@]}")"
    "$(printf '  %-20s %d' 'Incomplete:'     "$incomplete_genes")"
    "$SEP"
)

# ─── Print & optionally save report ───────────────────────────────────────────
printf '%s\n' "${report_lines[@]}"

if [[ -n "$REPORT_FILE" ]]; then
    printf '%s\n' "${report_lines[@]}" > "$REPORT_FILE"
    echo
    info "Report also written to: $REPORT_FILE"
fi
