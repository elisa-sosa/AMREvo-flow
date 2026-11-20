#!/usr/bin/env bash
# Materialise the pinned upstream GEMSWEEP checkout and compatibility patch.
set -euo pipefail

PROJECT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
CHECKOUT=${1:-"$PROJECT_DIR/vendor/gemsweep/checkout"}
EXPECTED_COMMIT=d06885f
UPSTREAM=https://github.com/sanger-pathogens/gemsweep.git
PATCH="$PROJECT_DIR/patches/gemsweep-nextflow25.patch"

if [[ ! -d "$CHECKOUT/.git" ]]; then
  git clone --recurse-submodules "$UPSTREAM" "$CHECKOUT"
fi

git -C "$CHECKOUT" fetch --tags origin
git -C "$CHECKOUT" checkout --detach "$EXPECTED_COMMIT"
git -C "$CHECKOUT" submodule update --init --recursive

actual=$(git -C "$CHECKOUT" rev-parse --short=7 HEAD)
[[ "$actual" == "$EXPECTED_COMMIT" ]] || {
  echo "ERROR: GEMSWEEP revision is $actual, expected $EXPECTED_COMMIT" >&2
  exit 2
}

if git -C "$CHECKOUT" apply --reverse --check "$PATCH" 2>/dev/null; then
  echo "GEMSWEEP compatibility patch is already applied."
else
  [[ -z "$(git -C "$CHECKOUT" status --short -- subworkflows/validate_prebuilt_input.nf)" ]] || {
    echo "ERROR: GEMSWEEP validation subworkflow has an unreviewed local change" >&2
    exit 2
  }
  git -C "$CHECKOUT" apply --check "$PATCH"
  git -C "$CHECKOUT" apply "$PATCH"
fi

git -C "$CHECKOUT" diff --check
test -s "$CHECKOUT/LICENSE.md"
printf 'Pinned GEMSWEEP checkout ready: %s (%s)\n' "$CHECKOUT" "$actual"
