#!/usr/bin/env bash


# Works in junction with the curtis2026_reproduction.py script

# Author: N.Musto
# Last updated: 17-10-2026


set -euo pipefail

# Run from the working folder containing the Python script.
project_dir="$PWD"
python_script="$project_dir/curtis2026_reproduction.py"
# First argument: seed. Remaining arguments: draw counts.
seed="${1:-2026}"

if (( $# > 1 )); then
    draw_counts=("${@:2}")
else
    draw_counts=(20 200 2000)
fi



if [[ ! -f "$python_script" ]]; then
    printf 'Cannot find: %s\nRun this Bash script from your project folder.\n' \
        "$python_script" >&2
    exit 1
fi

# Make a separate batch folder each time, even if runs start in the same second.
output_root="$project_dir/results/figureA3_batches"
mkdir -p -- "$output_root"
batch_dir=$(mktemp -d "$output_root/$(date +%Y%m%d_%H%M%S)_XXXXXX")

# Retain the exact code used for this batch, and run the Python snapshot.
cp -- "$python_script" "$batch_dir/curtis2026_reproduction.py"
cp -- "${BASH_SOURCE[0]}" "$batch_dir/test_runs_figa3.bash"

printf '\nBatch folder: %s\n' "$batch_dir"

for draws in "${draw_counts[@]}"; do
    run_name="figureA3_d${draws}_s${seed}"
    run_dir="$batch_dir/$run_name"
    mkdir -- "$run_dir"

    printf '\nRunning %s draws with seed %s...\n' "$draws" "$seed"

    # An array keeps paths containing spaces intact.
    command=(python3 "$batch_dir/curtis2026_reproduction.py"
        --figure A3
        --ratio 1.0008003201280513
        --noise-halfwidth 0.0002
        --max-initial 1000
        --draws "$draws"
        --show-draws 20
        --seed "$seed"
        --summary-only
        --output-prefix "$run_dir/$run_name")

    # Save the command, then show and record Python's output and errors.
    {
        printf 'Command: '
        printf '%q ' "${command[@]}"
        printf '\n\n'
    } > "$run_dir/run.log"

    "${command[@]}" 2>&1 | tee -a "$run_dir/run.log"
done

printf '\nAll runs completed. Results: %s\n' "$batch_dir"

