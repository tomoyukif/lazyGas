#!/usr/bin/env bash
# Build or pull the lazyGas Ollama Apptainer image.
set -euo pipefail

APPTAINER="${LAZYGAS_APPTAINER_BIN:-$(command -v apptainer 2>/dev/null || command -v singularity 2>/dev/null || true)}"
if [[ -z "${APPTAINER}" ]]; then
  echo "apptainer or singularity not found on PATH" >&2
  exit 1
fi

DEST="${1:-${LAZYGAS_OLLAMA_SIF:-$HOME/.cache/lazygas/ollama/ollama.sif}}"
METHOD="${LAZYGAS_OLLAMA_BUILD_METHOD:-pull}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEF_FILE="${SCRIPT_DIR}/ollama.def"

mkdir -p "$(dirname "${DEST}")"

if [[ "${METHOD}" == "build" ]]; then
  echo "Building ${DEST} from ${DEF_FILE} ..."
  "${APPTAINER}" build "${DEST}" "${DEF_FILE}"
else
  echo "Pulling docker://ollama/ollama:latest -> ${DEST} ..."
  "${APPTAINER}" pull "${DEST}" docker://ollama/ollama:latest
fi

echo "Done: ${DEST}"
