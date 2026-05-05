#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: ./run_omp_thread_sweep.sh <a_max> [thread_list]

Arguments:
  a_max        Search limit passed to ./bin/solve_omp
  thread_list  Thread counts to test, separated by ':' or ','
               Default: 1:2:4:8:16:32

Environment overrides:
  OMP_THREAD_LIST=1:2:4:8:16:32
  OMP_PROC_BIND=false
  OMP_PLACES=
USAGE
}

if [[ $# -lt 1 || $# -gt 2 ]]; then
  usage
  exit 2
fi

A_MAX="$1"
THREAD_LIST="${2:-${OMP_THREAD_LIST:-1:2:4:8:16:32}}"
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${PROJECT_DIR}/logs"
RESULTS_DIR="${PROJECT_DIR}/results"
STAMP="$(date +%Y%m%d_%H%M%S)"
CSV_PATH="${RESULTS_DIR}/omp_thread_sweep_a${A_MAX}_${STAMP}.csv"

mkdir -p "${LOG_DIR}" "${RESULTS_DIR}"

make -C "${PROJECT_DIR}" omp

printf 'threads,elapsed_ms,solution_count,status,log_file\n' > "${CSV_PATH}"

IFS=':,' read -r -a THREAD_COUNTS <<< "${THREAD_LIST}"
for threads in "${THREAD_COUNTS[@]}"; do
  LOG_PATH="${LOG_DIR}/omp-a${A_MAX}-t${threads}-${STAMP}.log"
  printf '=== OMP_NUM_THREADS=%s ===\n' "${threads}"

  set +e
  OMP_NUM_THREADS="${threads}" \
  OMP_PROC_BIND="${OMP_PROC_BIND:-false}" \
  ./bin/solve_omp "${A_MAX}" > "${LOG_PATH}" 2>&1
  status=$?
  set -e

  elapsed_ms="$(sed -n 's/.*Finished search in \([0-9][0-9]*\) ms!.*/\1/p' "${LOG_PATH}" | tail -n 1)"
  solution_count="$(grep -c '\^6 + .* = .*' "${LOG_PATH}" || true)"

  if [[ ${status} -eq 0 && -n "${elapsed_ms}" ]]; then
    printf '%s,%s,%s,ok,%s\n' \
      "${threads}" "${elapsed_ms}" "${solution_count}" "${LOG_PATH}" >> "${CSV_PATH}"
  else
    printf '%s,,%s,failed,%s\n' \
      "${threads}" "${solution_count}" "${LOG_PATH}" >> "${CSV_PATH}"
  fi
done

printf 'Saved results to %s\n' "${CSV_PATH}"
