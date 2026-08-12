#!/usr/bin/env bash
# scripts/branch_cleanup_audit.sh -- RYA-793
# Safe LOCAL-branch prune. Dry-run by default; on --execute deletes ONLY branches whose
# commits are provably in main (by ancestry, or by patch-equivalence for squash/rebase
# merges), logging each tip SHA first so every deletion is recoverable.
# NEVER touches: main, the current HEAD branch, the protected in-flight patterns, any
# unmerged branch, or any branch checked out in another worktree.
# Remote branches are NOT touched.
#
# ---------------------------------------------------------------------------------
# THREE CORRECTIONS TO THE ORIGINAL SPEC, each forced by measuring this repo
# ---------------------------------------------------------------------------------
# 1. `: gone]` FINDS NOTHING HERE. The spec detects squash-merges by looking for a
#    'gone' upstream. Measured on the Mac: 0 of 99 branches have one -- because the
#    remote still holds 246 heads, so no upstream has been deleted. GitHub's "auto
#    delete head branches" is off, which the spec's own Phase-2 note identifies as the
#    root cause. So the squash detector the script exists for would have matched zero
#    branches. Replaced with `git cherry` PATCH-EQUIVALENCE against main, which is the
#    RYA-568 ratified test and does not care whether the remote branch survives.
#    Effect: 3 squash/rebase-merged branches found that 'gone' could never have seen.
#
# 2. LOCAL `main` IS STALE, so it is the wrong authority (the RYA-510 gotcha).
#    Measured: `--merged main` = 30 branches, `--merged origin/main` = 51. Twenty-one
#    merged branches are invisible to the local ref, and `git branch -d` -- which
#    checks against HEAD, not against the remote -- would refuse them as "not fully
#    merged". This script classifies against origin/main after a fetch, and says so
#    loudly when the two refs disagree.
#
# 3. 36 OF 99 BRANCHES ARE CHECKED OUT IN WORKTREES. `git branch -d/-D` refuses those
#    outright, so listing them as DELETE candidates would promise deletions that cannot
#    happen. They get their own bucket. (Related: the spec's `git branch -vv | awk
#    '{print $1}'` would also mis-parse them -- git prefixes those lines with `+ `, and
#    the current branch with `* `, so $1 is the marker, not the branch name.)
# ---------------------------------------------------------------------------------
set -euo pipefail

MAIN_LOCAL="main"
MAIN="origin/main"                            # the authority; see correction 2
PROTECT_PATTERNS='rya-784|rya-785|rya-762'    # current in-flight work; extend as needed
LOG="branch_cleanup_$(date +%Y%m%d_%H%M%S).log"
EXECUTE=0
[[ "${1:-}" == "--execute" ]] && EXECUTE=1

git fetch --prune --quiet
CUR=$(git rev-parse --abbrev-ref HEAD)

if ! git rev-parse --verify --quiet "$MAIN" >/dev/null; then
  echo "FATAL: $MAIN does not resolve. Refusing to classify against a missing ref."
  exit 1
fi

# Branches checked out in ANY worktree (including this one). `git branch -d` cannot
# remove these, so they are reported, never offered for deletion.
# Newline-delimited string, not an associative array: macOS ships bash 3.2, which has
# no `declare -A`, and this script has to run on both machines.
WORKTREE_BRANCHES=$(git worktree list --porcelain \
                    | awk '/^branch /{sub("refs/heads/","",$2); print $2}')
in_worktree() {
  printf '%s\n' "$WORKTREE_BRANCHES" | grep -qxF "$1"
}

if [[ "$(git rev-parse --verify --quiet "$MAIN_LOCAL" || echo none)" \
      != "$(git rev-parse "$MAIN")" ]]; then
  echo "NOTE: local '$MAIN_LOCAL' != '$MAIN' -- classifying against $MAIN (correction 2)."
fi

is_protected() {
  local b="$1"
  [[ "$b" == "$MAIN_LOCAL" || "$b" == "$CUR" ]] && return 0
  [[ "$b" =~ $PROTECT_PATTERNS ]] && return 0
  return 1
}

declare -a DELETE=() KEEP=() REVIEW=() BLOCKED=() UNMERGED=()

while IFS= read -r b; do
  [[ -z "$b" ]] && continue

  # main falls THROUGH to is_protected rather than being skipped, so it appears in the
  # PROTECTED bucket. The spec's smoke test asks the reader to verify main is protected;
  # a branch that is silently never considered looks identical to one that was missed.
  if is_protected "$b"; then KEEP+=("$b"); continue; fi

  # Is every commit on this branch already in main?
  #   ANCESTOR      -- a plain merge; git's own -d will agree.
  #   IN-MAIN       -- not an ancestor, but `git cherry` reports no '+' commits, i.e.
  #                    every patch is upstream. That is what a squash or rebase merge
  #                    leaves behind, and it is the ONLY reliable way to see one.
  #   otherwise     -- real unmerged work. Never a candidate.
  if git merge-base --is-ancestor "$b" "$MAIN" 2>/dev/null; then
    why="merged"
  elif ! git cherry "$MAIN" "$b" 2>/dev/null | grep -q '^+'; then
    why="squash-merged"
  else
    # Unmerged. If its upstream is ALSO gone, that is the spec's review case: the
    # remote branch was deleted without the work landing. Surface it, never delete it.
    if [[ -n "$(git for-each-ref --format='%(upstream:track)' "refs/heads/$b")" ]] \
       && git for-each-ref --format='%(upstream:track)' "refs/heads/$b" | grep -q 'gone'; then
      REVIEW+=("$b|gone-upstream-but-commits-not-in-main")
    else
      UNMERGED+=("$b")
    fi
    continue
  fi

  if in_worktree "$b"; then
    BLOCKED+=("$b|$why|checked-out-in-a-worktree")
  else
    DELETE+=("$b|$why")
  fi
done < <(git for-each-ref --format='%(refname:short)' refs/heads/)

# bash 3.2 + `set -u`: expanding an empty array errors, so every bucket prints through
# this guard rather than relying on ${arr[@]:-default}.
print_bucket() {
  local title="$1"; shift
  printf '\n== %s (%d) ==\n' "$title" "$#"
  if [[ "$#" -eq 0 ]]; then printf '  <none>\n'; else printf '  %s\n' "$@"; fi
}

print_bucket "PROTECTED / KEPT" ${KEEP[@]+"${KEEP[@]}"}
print_bucket "UNMERGED, left alone" ${UNMERGED[@]+"${UNMERGED[@]}"}
print_bucket "MANUAL REVIEW -- gone upstream, commits NOT in main" \
             ${REVIEW[@]+"${REVIEW[@]}"}
print_bucket "MERGED but BLOCKED -- checked out in a worktree, cannot delete" \
             ${BLOCKED[@]+"${BLOCKED[@]}"}
print_bucket "DELETE CANDIDATES" ${DELETE[@]+"${DELETE[@]}"}

if [[ "$EXECUTE" -eq 0 ]]; then
  printf '\nDRY RUN -- nothing deleted. Review the candidates, then re-run with --execute.\n'
  exit 0
fi

printf '\nDeleting (tip SHAs -> %s for recovery)...\n' "$LOG"
for entry in ${DELETE[@]+"${DELETE[@]}"}; do
  [[ -z "$entry" ]] && continue
  b="${entry%%|*}"; why="${entry##*|}"
  sha=$(git rev-parse "$b")
  echo "$sha  $b  ($why)  -- recover: git branch $b $sha" | tee -a "$LOG"
  # -d first: git's own refuse-if-unmerged net. It checks against HEAD, so with a stale
  # local main it can refuse a branch we have ALREADY proven is an ancestor of
  # origin/main -- a strictly stronger test. Only then fall back to -D, and only for
  # branches that passed that proof above.
  if ! git branch -d "$b" 2>/dev/null; then
    git branch -D "$b"
    echo "    (-d refused -- deleted with -D; ancestry to $MAIN was proven above)" \
      | tee -a "$LOG"
  fi
done
printf 'Done. Recovery lines in %s (also in git reflog, ~90 days).\n' "$LOG"
