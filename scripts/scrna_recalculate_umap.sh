#!/usr/bin/env bash
#SBATCH --job-name=codespring_umap_view
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --export=NONE

set -euo pipefail

export PATH="${PATH:-/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"
export MKL_NUM_THREADS="$OMP_NUM_THREADS"
export NUMBA_NUM_THREADS="$OMP_NUM_THREADS"

if ! command -v module >/dev/null 2>&1; then
  for module_init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash; do
    if [[ -r "$module_init" ]]; then
      set +u
      # shellcheck disable=SC1090
      source "$module_init"
      set -u
      break
    fi
  done
fi

engine="${1:?engine is required}"
object_path="${2:?processed object is required}"
cells_path="${3:?included-cell table is required}"
output_path="${4:?output path is required}"
n_neighbors="${5:-15}"
min_dist="${6:-0.3}"
seed="${7:-1234}"
scanpy_container="${8:-}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$(dirname "$output_path")"
tmp_root="$(dirname "$output_path")/tmp"
mkdir -p "$tmp_root"
export TMPDIR="$tmp_root"
export TMP="$tmp_root"
export TEMP="$tmp_root"
export MPLCONFIGDIR="$tmp_root/matplotlib"
export NUMBA_CACHE_DIR="$tmp_root/numba"
export XDG_CACHE_HOME="$tmp_root/cache"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR" "$XDG_CACHE_HOME"

case "${engine,,}" in
  seurat)
    module purge >/dev/null 2>&1 || true
    module load EB5Modules >/dev/null 2>&1
    module load "${CSL_SEURAT_MODULE:-Seurat/5.4.0-foss-2024a-R-4.4.2}" >/dev/null 2>&1
    export R_LIBS_USER="$(dirname "$output_path")/.codespring_unused_user_library"
    export R_ENVIRON_USER=/dev/null
    export R_PROFILE_USER=/dev/null
    Rscript "$script_dir/scrna_recalculate_umap_seurat.R" "$object_path" "$cells_path" "$output_path" "$n_neighbors" "$min_dist" "$seed"
    ;;
  scanpy)
    module load EBModules >/dev/null 2>&1 || true
    module load singularity/3.6.3 >/dev/null 2>&1 || true
    if [[ -z "$scanpy_container" || ! -r "$scanpy_container" ]]; then
      echo "ERROR: The shared Scanpy container is unavailable: $scanpy_container" >&2
      exit 2
    fi
    singularity exec --cleanenv \
      --bind="$(dirname "$object_path")" \
      --bind="$(dirname "$cells_path")" \
      --bind="$script_dir" \
      "$scanpy_container" python "$script_dir/scrna_recalculate_umap_scanpy.py" \
      "$object_path" "$cells_path" "$output_path" "$n_neighbors" "$min_dist" "$seed"
    ;;
  *)
    echo "ERROR: engine must be seurat or scanpy" >&2
    exit 2
    ;;
esac

test -s "$output_path"
rm -rf "$tmp_root"
echo "Selected-cell UMAP ready: $output_path"
