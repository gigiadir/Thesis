#!/bin/bash
# Submit Stage 5 (gene-shuffle nulls) to SGE.
#
# Usage:
#   bash submit_stage05_nulls.sh
#   NCORES=16 bash submit_stage05_nulls.sh
#   SGE_PE=shared NCORES=16 bash submit_stage05_nulls.sh
#   SGE_PE=none bash submit_stage05_nulls.sh    # single-slot, no -pe
#
# Monitor:
#   qstat -u "$USER"
#
# Logs:
#   $OUTPUT_DIR/results/logs/atlas_nulls.o* / .e*

set -euo pipefail

ATLAS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HOME_DIR="${HOME:-/gpfs0/bgu-ofircohen/users/gigiadir}"
CONDA_ENV_PREFIX="${SCDIFFCOM_ENV_PREFIX:-${HOME_DIR}/.conda/envs/scDiffComPipeline_env}"
R_LIBS_SITE_EXPORT="${R_LIBS_SITE:-${CONDA_ENV_PREFIX}/lib/R/library}"

CONFIG_PATH="${CONFIG_PATH:-${ATLAS_DIR}/config_batch.yml}"
OUTPUT_DIR="${OUTPUT_DIR:-$HOME/CCC-scDiffCom/results/reproducibility_atlas/v1-v2}"

# Use same log root pattern as scDiffCom jobs (known writable on exec nodes)
LOG_DIR="${LOG_DIR:-$HOME/CCC-scDiffCom/results/reproducibility_atlas/logs}"
mkdir -p "${LOG_DIR}"

LOG_OUT="${LOG_DIR}/atlas_nulls.o"
LOG_ERR="${LOG_DIR}/atlas_nulls.e"

NCORES="${NCORES:-64}"
SGE_PE="${SGE_PE:-shared}"
N_PERM_ARG="${N_PERM_ARG:-}"
EXTRA_ARGS="${EXTRA_ARGS:-}"

RUN_SCRIPT="${ATLAS_DIR}/run_stage_05_nulls.R"
ENV_WRAPPER="${HOME_DIR}/Scripts/run_with_scDiffComPipeline_env.sh"

if [[ ! -x "${RUN_SCRIPT}" && ! -f "${RUN_SCRIPT}" ]]; then
  echo "ERROR: run_stage_05_nulls.R not found: ${RUN_SCRIPT}" >&2
  exit 1
fi

if [[ ! -f "${ENV_WRAPPER}" ]]; then
  echo "ERROR: env wrapper not found: ${ENV_WRAPPER}" >&2
  exit 1
fi

ARGS=(--config "${CONFIG_PATH}" --ncores "${NCORES}")
if [[ -n "${N_PERM_ARG}" ]]; then
  ARGS+=(--n-perm "${N_PERM_ARG}")
fi
if [[ -n "${EXTRA_ARGS}" ]]; then
  # shellcheck disable=SC2206
  ARGS+=(${EXTRA_ARGS})
fi

PE_ARGS=()
if [[ "${SGE_PE}" != "none" && -n "${SGE_PE}" ]]; then
  PE_ARGS=(-pe "${SGE_PE}" "${NCORES}")
fi

echo "Submitting atlas Stage 5 nulls"
echo "  queue:      intel_all.q"
echo "  parallel:   ${SGE_PE:-none}"
echo "  cores:      ${NCORES}"
echo "  atlas dir:  ${ATLAS_DIR}"
echo "  output dir: ${OUTPUT_DIR}"
echo "  logs:       ${LOG_DIR}/atlas_nulls.o* / .e*"

qsub -q intel_all.q -S /bin/bash -cwd \
  -N atlas_nulls \
  "${PE_ARGS[@]}" \
  -o "${LOG_OUT}" \
  -e "${LOG_ERR}" \
  -v "SCDIFFCOM_PIPELINE_R=${RUN_SCRIPT},OUTPUT_DIR=${OUTPUT_DIR},R_LIBS_SITE=${R_LIBS_SITE_EXPORT},SCDIFFCOM_ENV_PREFIX=${CONDA_ENV_PREFIX}" \
  "${ENV_WRAPPER}" \
  "${ARGS[@]}"
