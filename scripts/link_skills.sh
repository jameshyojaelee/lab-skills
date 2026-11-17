#!/usr/bin/env bash
# Link every skill in this repo into agent skill folders, so each tool reads this one copy.
#   scripts/link_skills.sh                      # default folders below
#   scripts/link_skills.sh ~/some/skills/dir    # or name folders
# Existing links are replaced; real folders with the same name are left alone and reported.
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
TARGETS=("$@")
[ ${#TARGETS[@]} -eq 0 ] && TARGETS=("$HOME/.claude/skills" "$HOME/.codex/skills" "$HOME/.agents/skills")
for t in "${TARGETS[@]}"; do
  mkdir -p "$t"
  for s in "$REPO"/skills/*/; do
    n="$(basename "$s")"
    if [ -e "$t/$n" ] && [ ! -L "$t/$n" ]; then
      echo "skip $t/$n (a real folder exists; move it away first)"; continue
    fi
    ln -sfn "${s%/}" "$t/$n" && echo "linked $t/$n"
  done
done
