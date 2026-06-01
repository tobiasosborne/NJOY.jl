#!/usr/bin/env bash
#
# setup_reference.sh — align the local Fortran reference oracle (njoy-reference/)
# to the SHA pinned in REFERENCE_PIN.
#
# njoy-reference/ is a git clone of NJOY2016 used as the canonical oracle for
# bit-identical validation. It is git-ignored by NJOY.jl and cloned per-machine.
# Upstream `develop` rebases, so we pin a fixed commit (see REFERENCE_PIN) to
# keep the oracle from moving under us.
#
# Idempotent: clones if missing, otherwise fetches and checks out the pin.
# Usage:
#     bash scripts/setup_reference.sh
set -euo pipefail

# Resolve repo root from this script's own location (works from any CWD).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

PIN_FILE="${REPO_ROOT}/REFERENCE_PIN"
REF_DIR="${REPO_ROOT}/njoy-reference"
URL="https://github.com/njoy/NJOY2016.git"

if [[ ! -f "${PIN_FILE}" ]]; then
    echo "FAIL: REFERENCE_PIN not found at ${PIN_FILE}" >&2
    exit 1
fi

# Parse the pinned SHA from the `SHA=...` line.
SHA="$(grep -E '^SHA=' "${PIN_FILE}" | head -n1 | cut -d= -f2 | tr -d '[:space:]')"
if [[ -z "${SHA}" ]]; then
    echo "FAIL: could not parse SHA= from ${PIN_FILE}" >&2
    exit 1
fi
echo "Pinned SHA: ${SHA}"

if [[ ! -d "${REF_DIR}" ]]; then
    echo "njoy-reference/ absent — cloning ${URL} ..."
    git clone "${URL}" "${REF_DIR}"
else
    echo "njoy-reference/ present — fetching origin ..."
    git -C "${REF_DIR}" fetch --quiet origin
fi

echo "Checking out ${SHA} ..."
if ! git -C "${REF_DIR}" checkout "${SHA}"; then
    echo "FAIL: could not checkout ${SHA} in ${REF_DIR}" >&2
    echo "      The pinned commit may not be present locally; try a fresh clone." >&2
    exit 1
fi

HEAD="$(git -C "${REF_DIR}" rev-parse HEAD)"
if [[ "${HEAD}" == "${SHA}" ]]; then
    echo "OK: njoy-reference HEAD is at the pin (${HEAD})"
else
    echo "FAIL: njoy-reference HEAD is ${HEAD}, expected ${SHA}" >&2
    exit 1
fi

echo
echo "Note: the Fortran oracle binary (njoy-reference/build/njoy) may need a"
echo "rebuild after this checkout if you use it for fresh oracle generation."
