#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_ROOT}/docs"

QUARTO_PYTHON_CANDIDATE="${REPO_ROOT}/.venv/bin/python"
if [[ ! -x "${QUARTO_PYTHON_CANDIDATE}" && -x "${REPO_ROOT}/../.venv/bin/python" ]]; then
  QUARTO_PYTHON_CANDIDATE="${REPO_ROOT}/../.venv/bin/python"
fi
export QUARTO_PYTHON="${QUARTO_PYTHON_CANDIDATE}"

rm -rf reference
mkdir -p reference

uv run quartodoc build --config _quarto.yml
quarto render
