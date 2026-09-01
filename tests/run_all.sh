#!/usr/bin/env bash
set -Eeuo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
lab_root="${1:-$(cd -- "$repo_root/../CodeSpringLab-fix" && pwd)}"
(cd "$repo_root" && bash tests/smoke_test_fetchngs.sh)
(cd "$repo_root" && Rscript --vanilla tests/test_sarek_manifest.R)
(cd "$repo_root" && Rscript --vanilla tests/test_sarek_bam_inspector.R)
(cd "$repo_root" && bash tests/test_sarek_samtools_launcher.sh)
(cd "$repo_root" && Rscript --vanilla tests/test_sarek_manifest_shiny.R)
(cd "$repo_root" && Rscript --vanilla tests/test_sarek_nextflow_input.R)
(cd "$repo_root" && Rscript --vanilla tests/test_sarek_sex_fallback.R)
(cd "$repo_root" && Rscript --vanilla tests/test_sarek_live_activity.R)
(cd "$repo_root" && Rscript --vanilla tests/test_sarek_submission.R)
(cd "$repo_root" && Rscript tests/smoke_test_app_helpers.R "$lab_root")
CSL_CODESPRINGLAB_ROOT="$lab_root" "$repo_root/run_codespringweb.sh" --check-config
echo "All CodeSpringApp tests passed."
