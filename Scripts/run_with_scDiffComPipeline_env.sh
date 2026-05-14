#!/usr/bin/env bash
# Run scDiffComPipeline.R inside the conda env scDiffComPipeline_env (batch-safe).
# Usage: same flags as Rscript would take after the pipeline script, e.g.
#   run_with_scDiffComPipeline_env.sh --gene ABI1 --dataset_path ... --patient_summary_path ...
#
# Optional overrides:
#   SCDIFFCOM_CONDA_ENV     conda env name (default: scDiffComPipeline_env)
#   SCDIFFCOM_ENV_PREFIX    full path to the env (default: $HOME/.conda/envs/$SCDIFFCOM_CONDA_ENV)
#   SCDIFFCOM_PIPELINE_R    absolute path to scDiffComPipeline.R
#   CONDA_SH                path to conda.sh — only used if $SCDIFFCOM_ENV_PREFIX/bin/Rscript is missing

set -euo pipefail

# Save CLI args before any function call: inside functions, "$@" is the *function* args, not the script's.
argv=("$@")

SCDIFFCOM_CONDA_ENV="${SCDIFFCOM_CONDA_ENV:-scDiffComPipeline_env}"
SCDIFFCOM_PIPELINE_R="${SCDIFFCOM_PIPELINE_R:-$HOME/Scripts/scDiffComPipeline.R}"
SCDIFFCOM_ENV_PREFIX="${SCDIFFCOM_ENV_PREFIX:-${HOME}/.conda/envs/${SCDIFFCOM_CONDA_ENV}}"

__fail_if_missing_pipeline() {
  local f="$1"
  if [[ ! -f "$f" ]]; then
    echo "ERROR ($(hostname)): pipeline script not visible here: $f" >&2
    echo "PWD=$(pwd) HOME=${HOME-}" >&2
    ls -la "$(dirname "$f")" 2>&1 | head -50 >&2 || true
    exit 1
  fi
}

# Prefer activating by env directory (no conda.sh): matches e.g.
#   /gpfs0/bgu-ofircohen/users/gigiadir/.conda/envs/scDiffComPipeline_env
__use_env_prefix() {
  local prefix="$1"
  if [[ ! -x "${prefix}/bin/Rscript" ]]; then
    return 1
  fi
  export CONDA_PREFIX="$(cd "${prefix}" && pwd)"
  export PATH="${CONDA_PREFIX}/bin:${PATH}"
  unset R_HOME R_LIBS R_LIBS_USER
  export R_LIBS_SITE="${CONDA_PREFIX}/lib/R/library"
  __fail_if_missing_pipeline "${SCDIFFCOM_PIPELINE_R}"
  exec Rscript "${SCDIFFCOM_PIPELINE_R}" "${argv[@]}"
}

__source_conda() {
  if [[ -n "${CONDA_SH:-}" && -f "$CONDA_SH" ]]; then
    # shellcheck source=/dev/null
    source "$CONDA_SH"
    return 0
  fi
  local c
  for c in \
    "${HOME}/miniconda3/etc/profile.d/conda.sh" \
    "${HOME}/mambaforge/etc/profile.d/conda.sh" \
    "${HOME}/anaconda3/etc/profile.d/conda.sh" \
    "${HOME}/micromamba/etc/profile.d/conda.sh"; do
    if [[ -f "$c" ]]; then
      # shellcheck source=/dev/null
      source "$c"
      return 0
    fi
  done
  if command -v conda >/dev/null 2>&1; then
    local base
    base="$(conda info --base 2>/dev/null || true)"
    if [[ -n "$base" && -f "${base}/etc/profile.d/conda.sh" ]]; then
      # shellcheck source=/dev/null
      source "${base}/etc/profile.d/conda.sh"
      return 0
    fi
  fi
  echo "ERROR: Could not use env at ${SCDIFFCOM_ENV_PREFIX} (no bin/Rscript) and could not find conda.sh." >&2
  echo "Set SCDIFFCOM_ENV_PREFIX to your env path, or set CONDA_SH to .../etc/profile.d/conda.sh" >&2
  return 1
}

if ! __use_env_prefix "${SCDIFFCOM_ENV_PREFIX}"; then
  __source_conda
  conda activate "${SCDIFFCOM_CONDA_ENV}"

  unset R_HOME R_LIBS R_LIBS_USER
  export R_LIBS_SITE="${CONDA_PREFIX}/lib/R/library"

  __fail_if_missing_pipeline "${SCDIFFCOM_PIPELINE_R}"
  exec Rscript "${SCDIFFCOM_PIPELINE_R}" "${argv[@]}"
fi
