#!/usr/bin/env bash
set -Eeuo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
test_root="$(mktemp -d "${TMPDIR:-/tmp}/sarek-samtools-launcher.XXXXXX")"
trap 'rm -rf -- "$test_root"' EXIT

fake_bin="$test_root/bin"
module_log="$test_root/module.log"
module_init="$test_root/module-init.sh"
mkdir -p "$fake_bin"

cat > "$fake_bin/samtools" <<'SAMTOOLS'
#!/usr/bin/env bash
printf 'fake-samtools:%s\n' "$*"
SAMTOOLS
chmod 0755 "$fake_bin/samtools"

cat > "$module_init" <<'MODULE_INIT'
module() {
  if [[ "$1" != "load" ]]; then
    return 2
  fi
  printf '%s\n' "$2" >> "$CSL_SAREK_TEST_MODULE_LOG"
  case "$2" in
    EB5)
      return 0
      ;;
    SAMtools/1.23.1-GCC-15.2.0)
      PATH="$CSL_SAREK_TEST_FAKE_BIN:$PATH"
      export PATH
      return 0
      ;;
    *)
      return 1
      ;;
  esac
}
MODULE_INIT

output="$({
  CSL_SAREK_MODULE_INIT="$module_init" \
  CSL_SAREK_TEST_MODULE_LOG="$module_log" \
  CSL_SAREK_TEST_FAKE_BIN="$fake_bin" \
    "$repo_root/bin/sarek-samtools" quickcheck -v /tmp/example.bam
})"

[[ "$output" == "fake-samtools:quickcheck -v /tmp/example.bam" ]]
mapfile -t loaded_modules < "$module_log"
[[ "${loaded_modules[0]}" == "EB5" ]]
[[ "${loaded_modules[1]}" == "SAMtools/1.23.1-GCC-15.2.0" ]]
[[ "${#loaded_modules[@]}" -eq 2 ]]

echo "Sarek samtools module launcher tests: PASS"
