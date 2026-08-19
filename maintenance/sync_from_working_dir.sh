#!/usr/bin/env bash
# Re-sync this publication repository from the internal analysis working
# directory.
#
# Why this is not a plain `git pull`: this repo's history was rewritten with
# git-filter-repo to drop a 451 MB intermediate matrix that exceeds GitHub's
# file-size limit, so its commit SHAs differ from the working directory's.
#
# git-filter-repo is deterministic: filtering the same commits with the same
# arguments always yields the same SHAs (verified 2026-08-19). So re-filtering
# an updated working directory reproduces this repo's history exactly, plus
# whatever commits are new — which makes the packaging commits replayable with
# a normal rebase.
#
# Usage:  maintenance/sync_from_working_dir.sh [/path/to/PeakATail_wd]
# Then review the result and push yourself. This script never pushes.

set -euo pipefail

WD="${1:-/mnt/ssd1/Projects/PeakATail_wd}"
PUB="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# Paths stripped from all history. Both spellings are needed: the files lived
# at amirtest/... before the 2026-08-12 reorganisation moved them under
# archive/. tools/PeakATail was a gitlink with no .gitmodules.
FILTER_ARGS=(--invert-paths --path archive/amirtest --path amirtest --path tools/PeakATail)

command -v git-filter-repo >/dev/null || {
  echo "git-filter-repo not found. Install with: pip install --user git-filter-repo" >&2
  exit 1
}

echo "==> Cloning working directory from $WD"
git clone --no-local --no-checkout "$WD" "$WORK/wd" -q

echo "==> Filtering oversized history"
git -C "$WORK/wd" filter-repo "${FILTER_ARGS[@]}" --force >/dev/null

UPSTREAM_TIP="$(git -C "$WORK/wd" rev-parse reorg-manuscript 2>/dev/null \
             || git -C "$WORK/wd" rev-parse HEAD)"
echo "==> Filtered upstream tip: $UPSTREAM_TIP"

echo "==> Fetching it into this repo"
git -C "$PUB" fetch -q "$WORK/wd" "$UPSTREAM_TIP":refs/heads/_synced_upstream --force

# Everything on main that is NOT in the filtered upstream is packaging work.
BASE="$(git -C "$PUB" merge-base main _synced_upstream)"
echo "==> Packaging commits to replay:"
git -C "$PUB" log --oneline "$BASE"..main | sed 's/^/      /'

echo "==> Rebasing packaging commits onto the new upstream"
git -C "$PUB" rebase --onto _synced_upstream "$BASE" main

git -C "$PUB" branch -D _synced_upstream -q

cat <<EOF

==> Done. History rebuilt on this branch.

    Still to refresh by hand if the analysis moved:
      source_data/   cp \$WD/results/figures/manuscript/*.tsv source_data/
      envs/          conda env export -n bench_<tool> --no-builds > envs/bench_<tool>.yml

    Review, then push:
      git -C "$PUB" log --oneline | head
      git -C "$PUB" push origin main
EOF
