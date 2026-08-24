#!/usr/bin/env bash

set -Eeuo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
test_root="$(mktemp -d /tmp/codespring-fetchngs-smoke.XXXXXX)"
test_home="$test_root/home"
input_file="$test_root/accessions.txt"

cleanup() {
  [[ "$test_root" == /tmp/codespring-fetchngs-smoke.* ]] && rm -rf -- "$test_root"
}
trap cleanup EXIT

mkdir -p "$test_home"
printf 'SRR14593545\nERR1160846\n' > "$input_file"

HOME="$test_home" "$repo_root/bin/codespringflow" fetchngs \
  --input "$input_file" \
  --name smoke_fastq \
  --dry-run

fastq_run="$test_home/csl_results/fetchngs/smoke_fastq"
test -d "$fastq_run/results"
test -d "$test_home/.codespringflow/cache/singularity"
test -d "$test_home/.codespringflow/work/fetchngs/smoke_fastq"
test ! -e "$fastq_run/job_id.txt"
test -f "$fastq_run/input/accessions.csv"
cmp "$input_file" "$fastq_run/input/accessions.csv"
grep -Fq "input: '$fastq_run/input/accessions.csv'" "$fastq_run/params.yml"
bash -n "$fastq_run/run.sbatch"
grep -Fq 'resume_args=(-name "smoke_fastq")' "$fastq_run/run.sbatch"
grep -Fq '    resume_args=(-resume)' "$fastq_run/run.sbatch"
! grep -Fq '    -name "smoke_fastq"' "$fastq_run/run.sbatch"
grep -Fq "download_method: 'sratools'" "$fastq_run/params.yml"
grep -Fq "ena_metadata_fields: 'run_accession," "$fastq_run/params.yml"
! grep -Fq "parent_study" "$fastq_run/params.yml"
grep -Fq "skip_fastq_download: false" "$fastq_run/params.yml"
grep -Fq $'fetchngs_version\t1.12.0' "$fastq_run/run_manifest.tsv"
grep -Fq $'ena_metadata_fields\trun_accession,' "$fastq_run/run_manifest.tsv"
grep -Fq $'nextflow_version\t24.04.4' "$fastq_run/run_manifest.tsv"

HOME="$test_home" "$repo_root/bin/codespringflow" fetchngs \
  --input "$input_file" \
  --name smoke_metadata \
  --metadata-only \
  --dry-run

metadata_run="$test_home/csl_results/fetchngs/smoke_metadata"
grep -Fq "skip_fastq_download: true" "$metadata_run/params.yml"
grep -Fq $'metadata_only\ttrue' "$metadata_run/run_manifest.tsv"

custom_results_root="$test_root/custom-fetchngs-results"
CODESPRINGFLOW_RESULTS_ROOT="$custom_results_root" HOME="$test_home" \
  "$repo_root/bin/codespringflow" fetchngs \
  --input "$input_file" \
  --name smoke_custom_root \
  --dry-run
custom_run="$custom_results_root/smoke_custom_root"
test -f "$custom_run/input/accessions.csv"
grep -Fq "outdir: '$custom_run/results'" "$custom_run/params.yml"
test "$(CODESPRINGFLOW_RESULTS_ROOT="$custom_results_root" HOME="$test_home" "$repo_root/bin/codespringflow" root)" = "$custom_results_root"
custom_work_dir="$(awk -F $'\t' '$1 == "work_directory" { print $2 }' "$custom_run/run_manifest.tsv")"
custom_work_namespace="$(awk -F $'\t' '$1 == "work_namespace" { print $2 }' "$custom_run/run_manifest.tsv")"
[[ "$custom_work_namespace" == root_* ]]
[[ "$custom_work_dir" == "$test_home/.codespringflow/work/fetchngs/$custom_work_namespace/smoke_custom_root" ]]
test -d "$custom_work_dir"
test "$custom_work_dir" != "$test_home/.codespringflow/work/fetchngs/smoke_custom_root"

if HOME="$test_home" "$repo_root/bin/codespringflow" fetchngs --input "$input_file" --name smoke_fastq --dry-run >/dev/null 2>&1; then
  echo "Duplicate FetchNGS run names were not rejected." >&2
  exit 1
fi

if HOME="$test_home" "$repo_root/bin/codespringflow" fetchngs --input "$input_file" --name '../unsafe' --dry-run >/dev/null 2>&1; then
  echo "Unsafe FetchNGS run names were not rejected." >&2
  exit 1
fi

test "$(HOME="$test_home" "$repo_root/bin/codespringflow" root)" = "$test_home/csl_results/fetchngs"
test "$(HOME="$test_home" "$repo_root/bin/codespringflow" runtime-root)" = "$test_home/.codespringflow"

fake_bin="$test_root/fake-bin"
mkdir -p "$fake_bin"
printf '#!/usr/bin/env bash\nprintf "Submitted batch job 424242\\n"\n' > "$fake_bin/sbatch"
chmod 0755 "$fake_bin/sbatch"
printf '#!/usr/bin/env bash\nexit 0\n' > "$fake_bin/squeue"
printf '#!/usr/bin/env bash\nexit 0\n' > "$fake_bin/sacct"
chmod 0755 "$fake_bin/squeue" "$fake_bin/sacct"

# Execute a generated controller script with fake cluster/Nextflow commands and
# verify that initial and resumed runs receive mutually exclusive CLI options.
printf '#!/usr/bin/env bash\nexit 0\n' > "$fake_bin/module"
printf '#!/usr/bin/env bash\nprintf "%%s\\n" "$@" > "$NEXTFLOW_ARGS_FILE"\n' > "$fake_bin/nextflow"
chmod 0755 "$fake_bin/module" "$fake_bin/nextflow"
PATH="$fake_bin:$PATH" HOME="$test_home" NEXTFLOW_BIN="$fake_bin/nextflow" \
  "$repo_root/bin/codespringflow" fetchngs \
  --input "$input_file" \
  --name smoke_cli \
  --dry-run
cli_run="$test_home/csl_results/fetchngs/smoke_cli"
initial_args="$test_root/nextflow-initial.args"
resume_args="$test_root/nextflow-resume.args"
PATH="$fake_bin:$PATH" NEXTFLOW_ARGS_FILE="$initial_args" SLURM_JOB_ID=1001 \
  bash "$cli_run/run.sbatch"
grep -Fxq -- '-name' "$initial_args"
grep -Fxq -- 'smoke_cli' "$initial_args"
! grep -Fxq -- '-resume' "$initial_args"
PATH="$fake_bin:$PATH" NEXTFLOW_ARGS_FILE="$resume_args" SLURM_JOB_ID=1002 \
  CODESPRINGFLOW_RESUME=true bash "$cli_run/run.sbatch"
grep -Fxq -- '-resume' "$resume_args"
! grep -Fxq -- '-name' "$resume_args"

PATH="$fake_bin:$PATH" HOME="$test_home" NEXTFLOW_BIN=/bin/true \
  "$repo_root/bin/codespringflow" fetchngs \
  --input "$input_file" \
  --name smoke_submit

submitted_run="$test_home/csl_results/fetchngs/smoke_submit"
test "$(<"$submitted_run/job_id.txt")" = "424242"
test "$(wc -l < "$submitted_run/job_history.txt")" -eq 1

# Simulate a bundle generated before the ENA and resume-CLI compatibility fixes.
sed -i '/^ena_metadata_fields:/d' "$submitted_run/params.yml"
sed -i 's/^resume_args=(-name "smoke_submit")$/resume_args=()/' "$submitted_run/run.sbatch"
sed -i '/    -with-report/i\    -name "smoke_submit" \\' "$submitted_run/run.sbatch"
touch -d '10 minutes ago' "$submitted_run/job_id.txt"
PATH="$fake_bin:$PATH" HOME="$test_home" NEXTFLOW_BIN=/bin/true \
  "$repo_root/bin/codespringflow" resume smoke_submit
test "$(wc -l < "$submitted_run/job_history.txt")" -eq 2
grep -Fq "ena_metadata_fields: 'run_accession," "$submitted_run/params.yml"
! grep -Fq "parent_study" "$submitted_run/params.yml"
test "$(grep -c '^ena_metadata_fields:' "$submitted_run/params.yml")" -eq 1
test -f "$submitted_run/run.sbatch.pre-resume-cli-fix"
grep -Fq 'resume_args=(-name "smoke_submit")' "$submitted_run/run.sbatch"
! grep -Fq '    -name "smoke_submit"' "$submitted_run/run.sbatch"

# Reproduce the real legacy layout where the fixed -name option was embedded
# in the same physical line as the rest of the Nextflow command.
HOME="$test_home" "$repo_root/bin/codespringflow" fetchngs \
  --input "$input_file" \
  --name smoke_inline_legacy \
  --dry-run
inline_legacy_run="$test_home/csl_results/fetchngs/smoke_inline_legacy"
sed -i 's/^resume_args=(-name "smoke_inline_legacy")$/resume_args=()/' "$inline_legacy_run/run.sbatch"
sed -i 's/    -with-report/    -name "smoke_inline_legacy"    -with-report/' "$inline_legacy_run/run.sbatch"
PATH="$fake_bin:$PATH" HOME="$test_home" NEXTFLOW_BIN=/bin/true \
  "$repo_root/bin/codespringflow" resume smoke_inline_legacy
bash -n "$inline_legacy_run/run.sbatch"
grep -Fq 'resume_args=(-name "smoke_inline_legacy")' "$inline_legacy_run/run.sbatch"
! grep -Fq '    -name "smoke_inline_legacy"' "$inline_legacy_run/run.sbatch"
grep -Fq -- '-with-report' "$inline_legacy_run/run.sbatch"

# Recover an empty active script only from a validated non-empty backup.
HOME="$test_home" "$repo_root/bin/codespringflow" fetchngs \
  --input "$input_file" \
  --name smoke_empty_recovery \
  --dry-run
empty_recovery_run="$test_home/csl_results/fetchngs/smoke_empty_recovery"
cp -p -- "$empty_recovery_run/run.sbatch" "$empty_recovery_run/run.sbatch.pre-resume-cli-fix"
: > "$empty_recovery_run/run.sbatch"
PATH="$fake_bin:$PATH" HOME="$test_home" NEXTFLOW_BIN=/bin/true \
  "$repo_root/bin/codespringflow" resume smoke_empty_recovery
test -s "$empty_recovery_run/run.sbatch"
bash -n "$empty_recovery_run/run.sbatch"
grep -Fq 'resume_args=(-name "smoke_empty_recovery")' "$empty_recovery_run/run.sbatch"

# An empty script without a usable backup must never reach sbatch.
HOME="$test_home" "$repo_root/bin/codespringflow" fetchngs \
  --input "$input_file" \
  --name smoke_empty_reject \
  --dry-run
empty_reject_run="$test_home/csl_results/fetchngs/smoke_empty_reject"
: > "$empty_reject_run/run.sbatch"
empty_reject_error="$test_root/empty-reject.err"
if PATH="$fake_bin:$PATH" HOME="$test_home" NEXTFLOW_BIN=/bin/true \
  "$repo_root/bin/codespringflow" resume smoke_empty_reject >"$empty_reject_error" 2>&1; then
  echo "An empty FetchNGS run script without a backup was submitted." >&2
  exit 1
fi
grep -Fq "run script is empty and no non-empty recovery backup exists" "$empty_reject_error"
test ! -e "$empty_reject_run/job_id.txt"

# A scheduler-active job must block a second resume submission.
printf '#!/usr/bin/env bash\nprintf "RUNNING\\n"\n' > "$fake_bin/squeue"
printf '#!/usr/bin/env bash\nexit 0\n' > "$fake_bin/sacct"
chmod 0755 "$fake_bin/squeue" "$fake_bin/sacct"
active_history_count="$(wc -l < "$submitted_run/job_history.txt")"
active_error="$test_root/active-resume.err"
if PATH="$fake_bin:$PATH" HOME="$test_home" NEXTFLOW_BIN=/bin/true \
  "$repo_root/bin/codespringflow" resume smoke_submit >"$active_error" 2>&1; then
  echo "An active FetchNGS run accepted a duplicate resume submission." >&2
  exit 1
fi
grep -Fq "still has active SLURM job 424242 (RUNNING)" "$active_error"
test "$(wc -l < "$submitted_run/job_history.txt")" -eq "$active_history_count"

# If Slurm has not registered a fresh job yet, its recent job-id file still
# blocks a rapid second click during the scheduler visibility delay.
printf '#!/usr/bin/env bash\nexit 0\n' > "$fake_bin/squeue"
recent_error="$test_root/recent-resume.err"
if PATH="$fake_bin:$PATH" HOME="$test_home" NEXTFLOW_BIN=/bin/true \
  "$repo_root/bin/codespringflow" resume smoke_submit >"$recent_error" 2>&1; then
  echo "A recently submitted FetchNGS run accepted a duplicate resume." >&2
  exit 1
fi
grep -Fq "was submitted recently as SLURM job 424242" "$recent_error"
test "$(wc -l < "$submitted_run/job_history.txt")" -eq "$active_history_count"

# A successfully completed run must reject Resume before reaching sbatch.
printf '#!/usr/bin/env bash\nexit 0\n' > "$fake_bin/squeue"
printf '#!/usr/bin/env bash\nprintf "COMPLETED\\n"\n' > "$fake_bin/sacct"
chmod 0755 "$fake_bin/squeue" "$fake_bin/sacct"

completed_history_count="$(wc -l < "$submitted_run/job_history.txt")"
completed_error="$test_root/completed-resume.err"

if PATH="$fake_bin:$PATH" HOME="$test_home" NEXTFLOW_BIN=/bin/true \
  "$repo_root/bin/codespringflow" resume smoke_submit >"$completed_error" 2>&1; then
  echo "A completed FetchNGS run incorrectly accepted Resume." >&2
  exit 1
fi

grep -Fq "completed successfully as SLURM job 424242" "$completed_error"
grep -Fq "Resume is only available for failed or stopped runs" "$completed_error"
test "$(wc -l < "$submitted_run/job_history.txt")" -eq "$completed_history_count"

# A failed terminal state must remain resumable.
printf '#!/usr/bin/env bash\nprintf "FAILED\\n"\n' > "$fake_bin/sacct"
chmod 0755 "$fake_bin/sacct"

PATH="$fake_bin:$PATH" HOME="$test_home" NEXTFLOW_BIN=/bin/true \
  "$repo_root/bin/codespringflow" resume smoke_submit >"$test_root/failed-resume.out"

grep -Fq "Submitted batch job 424242" "$test_root/failed-resume.out"
test "$(wc -l < "$submitted_run/job_history.txt")" -eq "$((completed_history_count + 1))"

echo "CodeSpringApp FetchNGS bundle smoke tests passed."
