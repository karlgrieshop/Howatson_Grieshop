#!/bin/bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: bash fly_human_conservation/setup_environment.sh [--python /path/to/python3.12]

Creates fly_human_conservation/.venv_local and installs the pinned requirements.
Uses python3.12 on PATH unless --python is provided. Python itself is not installed.
Refuses to overwrite an existing environment.
EOF
}

python_bin="python3.12"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --python)
            if [[ $# -lt 2 || -z "$2" || "$2" == --* ]]; then
                echo "Error: --python requires a Python 3.12 executable." >&2
                exit 2
            fi
            python_bin="$2"
            shift 2
            ;;
        --help|-h) usage; exit 0 ;;
        *) echo "Error: unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

cd "$(dirname "$0")/.."
if ! command -v "$python_bin" >/dev/null 2>&1; then
    echo "Error: Python 3.12 not found: $python_bin. Install it or pass --python /path/to/python3.12." >&2
    exit 1
fi
if ! "$python_bin" -c 'import sys; sys.exit(sys.version_info[:2] != (3, 12))'; then
    echo "Error: $python_bin is not Python 3.12." >&2
    exit 1
fi
environment="fly_human_conservation/.venv_local"
if [[ -e "$environment" ]]; then
    echo "Error: $environment already exists; refusing to overwrite it." >&2
    exit 1
fi
"$python_bin" -m venv "$environment"
"$environment/bin/python" -m pip --isolated --disable-pip-version-check install \
    --only-binary=:all: --requirement fly_human_conservation/alignment_requirements.txt \
    --report "$environment/install_report.json"
"$environment/bin/python" -m pip --isolated check
echo "Local environment ready: $environment"
