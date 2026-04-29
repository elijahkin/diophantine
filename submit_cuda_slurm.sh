#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: ./submit_cuda_slurm.sh <a_max> [max_solutions]

Environment overrides:
  CUDA_MODULE=11.8.0
  CUDA_ARCH=sm_70
  SLURM_ACCOUNT=<account>
  SLURM_PARTITION=<partition>
  SLURM_QOS=<qos>
  SLURM_TIME=01:00:00
  SLURM_MEM=16G
  SLURM_CPUS_PER_TASK=4
  SLURM_GPUS=1
  SLURM_GRES=gpu:1
USAGE
}

if [[ $# -lt 1 || $# -gt 2 ]]; then
  usage
  exit 2
fi

A_MAX="$1"
MAX_SOLUTIONS="${2:-65536}"
CUDA_MODULE="${CUDA_MODULE:-11.8.0}"
CUDA_ARCH="${CUDA_ARCH:-sm_70}"
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "${PROJECT_DIR}/logs"

SBATCH_ARGS=(
  --job-name=diophantine-cuda
  --output="${PROJECT_DIR}/logs/diophantine-cuda-%j.out"
  --error="${PROJECT_DIR}/logs/diophantine-cuda-%j.err"
  --time="${SLURM_TIME:-01:00:00}"
  --mem="${SLURM_MEM:-16G}"
  --cpus-per-task="${SLURM_CPUS_PER_TASK:-4}"
)

if [[ -n "${SLURM_ACCOUNT:-}" ]]; then
  SBATCH_ARGS+=(--account="${SLURM_ACCOUNT}")
fi
if [[ -n "${SLURM_PARTITION:-}" ]]; then
  SBATCH_ARGS+=(--partition="${SLURM_PARTITION}")
fi
if [[ -n "${SLURM_QOS:-}" ]]; then
  SBATCH_ARGS+=(--qos="${SLURM_QOS}")
fi
if [[ -n "${SLURM_GRES:-}" ]]; then
  SBATCH_ARGS+=(--gres="${SLURM_GRES}")
else
  SBATCH_ARGS+=(--gres="gpu:${SLURM_GPUS:-1}")
fi

EXPORTS="ALL,A_MAX=${A_MAX},MAX_SOLUTIONS=${MAX_SOLUTIONS},PROJECT_DIR=${PROJECT_DIR},CUDA_MODULE=${CUDA_MODULE},CUDA_ARCH=${CUDA_ARCH}"
SBATCH_ARGS+=(--export="${EXPORTS}")

sbatch "${SBATCH_ARGS[@]}" <<'SBATCH'
#!/usr/bin/env bash
set -euo pipefail

cd "${PROJECT_DIR}"

if ! command -v module >/dev/null 2>&1; then
  # Some Slurm environments do not initialize environment modules for
  # non-interactive shells.
  source /etc/profile.d/modules.sh 2>/dev/null || true
fi

module purge
module load "cuda/${CUDA_MODULE}"

make cuda CUDA_ARCH="${CUDA_ARCH}" NVCC="$(command -v nvcc)"
./bin/solve_cuda "${A_MAX}" "${MAX_SOLUTIONS}"
SBATCH
