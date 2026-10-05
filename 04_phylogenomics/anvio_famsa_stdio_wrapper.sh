#!/bin/bash
# Title          : anvio_famsa_stdio_wrapper.sh
# Description    : FAMSA wrapper used by anvi-get-sequences-for-hmm-hits (called by 03_build_bacteria71_phylogenomic_tree.sh)
# Author         : Julius Eyiuche Nweze
# Date           : 2026/04/23
# Usage          : ./anvio_famsa_stdio_wrapper.sh <input.fa> <output.fa>

set -euo pipefail

REAL_FAMSA="${REAL_FAMSA:-$HOME/miniconda3/envs/anvio-8/bin/famsa}"

if [[ ! -x "$REAL_FAMSA" ]]; then
    echo "Missing executable famsa binary: $REAL_FAMSA" >&2
    exit 1
fi

if [[ $# -lt 2 ]]; then
    exec "$REAL_FAMSA" "$@"
fi

args=("$@")
input_arg="${args[0]}"
output_arg="${args[1]}"
remaining_args=("${args[@]:2}")
tmp_input=""
tmp_output=""
tmp_stderr=""

cleanup() {
    [[ -n "$tmp_input" && -f "$tmp_input" ]] && rm -f "$tmp_input"
    [[ -n "$tmp_output" && -f "$tmp_output" ]] && rm -f "$tmp_output"
    [[ -n "$tmp_stderr" && -f "$tmp_stderr" ]] && rm -f "$tmp_stderr"
}

trap cleanup EXIT

if [[ "$input_arg" == "STDIN" ]]; then
    tmp_input="$(mktemp)"
    cat > "$tmp_input"
    input_arg="$tmp_input"
fi

if [[ "$output_arg" == "STDOUT" ]]; then
    tmp_output="$(mktemp)"
    output_arg="$tmp_output"
fi

tmp_stderr="$(mktemp)"

if ! "$REAL_FAMSA" "$input_arg" "$output_arg" "${remaining_args[@]}" 2>"$tmp_stderr"; then
    cat "$tmp_stderr" >&2
    exit 1
fi

if [[ -n "$tmp_output" ]]; then
    printf 'FAMSA\n\n'
    cat "$tmp_output"
    printf '\nDone!\n'
fi
