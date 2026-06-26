#!/usr/bin/env bash
set -euo pipefail

UPSTREAM_REPO="${UPSTREAM_REPO:-https://github.com/joelnitta/SORTER2.git}"
UPSTREAM_REF="${UPSTREAM_REF:-rpackage}"
PKG_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY_DST="${PKG_ROOT}/inst/python"

CORE_SCRIPTS=(
  SORTER2_FormatReads.py
  SORTER2_Stage1A_TrimSPAdes.py
  SORTER2_Stage1B_AssembleOrthologs.py
  SORTER2_Stage2_PhaseOrthologs.py
  SORTER2_Stage3_PhaseHybrids.py
  SORTER2_Processor.py
  SORTER2_ProgenitorProcessor.py
)

TMP_DIR="$(mktemp -d)"
cleanup() {
  rm -rf "${TMP_DIR}"
}
trap cleanup EXIT

echo "Cloning ${UPSTREAM_REPO} (${UPSTREAM_REF})"
git clone --depth 1 --branch "${UPSTREAM_REF}" "${UPSTREAM_REPO}" "${TMP_DIR}"

mkdir -p "${PY_DST}"

for script in "${CORE_SCRIPTS[@]}"; do
  src="${TMP_DIR}/${script}"
  if [[ ! -f "${src}" ]]; then
    echo "Missing upstream script: ${script}" >&2
    exit 1
  fi
  cp "${src}" "${PY_DST}/${script}"
  echo "Synced ${script}"
done

upstream_sha="$(git -C "${TMP_DIR}" rev-parse HEAD)"
echo "Upstream commit: ${upstream_sha}"
echo "Sync complete. Review with: git -C "${PKG_ROOT}" status --short"
