#!/usr/bin/env bash
# Reproduce every numerical claim of
# "On the function field analogue of the Goormaghtigh equation".
# Requires `uv` (https://docs.astral.sh/uv/); deps resolve from each script's
# inline PEP 723 block. For a plain env instead: pip install -r requirements.txt
# then replace `uv run` with `python`.
set -euo pipefail
cd "$(dirname "$0")"

run() { echo; echo "==== $1 ===="; uv run "src/$1"; }

run genus_verification.py
run monodromy_verification.py
run disjointness_resultant.py
run enestrom_kakeya_check.py

echo
echo "All verifications finished."
