#!/bin/bash

set -m   # enable job control

cleanup() {
    echo "🛑 Batch interrupted, killing all child processes"
    kill -- -$$ 2>/dev/null
    exit 130
}
trap cleanup INT TERM

if [ $# -lt 1 ]; then
    echo "Usage: $0 <folder> [extra loop_weaver options, e.g. -no_dynamic -thinplate]"
    exit 1
fi

FOLDER="$1"
shift
EXTRA_OPTS=("$@")

if [ ! -d "$FOLDER" ]; then
    echo "Error: '$FOLDER' is not a directory"
    exit 1
fi

# Locate the loop_weaver executable relative to this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="$SCRIPT_DIR/../build"
WEAVER="$BUILD_DIR/loop_weaver"

if [ ! -x "$WEAVER" ]; then
    echo "Error: loop_weaver executable not found at $WEAVER (build it first)"
    exit 1
fi

# Resolve the input folder to an absolute path, then work from build:
# loop_weaver writes its outputs to the current directory
FOLDER="$(cd "$FOLDER" && pwd)"
cd "$BUILD_DIR"

TIMEOUT=1800
TIMEOUT_LOG="timeouts.txt"
CALL_DIR="$(pwd)"
DONE_DIR="$CALL_DIR/done"

# Options of this batch. The other parameters come from
# loopweaver_settings.txt next to loop_weaver (or its compiled defaults);
# extra options given after the folder are appended and win.
OPTS=(-no_symm -angle 45 -error_setup 1 -subsampl 2)

mkdir -p "$DONE_DIR"

echo "=== Run started at $(date): ${OPTS[*]} ${EXTRA_OPTS[*]} ===" >> "$TIMEOUT_LOG"

shopt -s nullglob
for file in "$FOLDER"/*.obj; do
    echo "Processing: $file"

    start_time=$(date +%s)
    gtimeout --kill-after=10s "$TIMEOUT" \
        "$WEAVER" "$file" -batch "${OPTS[@]}" "${EXTRA_OPTS[@]}"

    status=$?
    end_time=$(date +%s)
    elapsed=$((end_time - start_time))

    if [ $status -eq 124 ]; then
        echo "⏱️  Timed out: $file (${elapsed}s)"
        echo "$(date '+%Y-%m-%d %H:%M:%S')  TIMEOUT  $file  (${elapsed}s)" >> "$TIMEOUT_LOG"
        continue
    elif [ $status -ne 0 ]; then
        echo "❌ Error ($status): $file (${elapsed}s)"
        echo "$(date '+%Y-%m-%d %H:%M:%S')  ERROR    $file  (${elapsed}s)" >> "$TIMEOUT_LOG"
        continue
    fi

    # Move all output files whose name contains the input stem
    stem=$(basename "$file" .obj)
    moved=0
    for out in "$CALL_DIR"/*"${stem}"*; do
        # Skip the done dir itself and the timeout log
        [ -d "$out" ] && continue
        [ "$out" = "$CALL_DIR/$TIMEOUT_LOG" ] && continue
        mv "$out" "$DONE_DIR/"
        moved=$((moved + 1))
    done

    echo "✅ Success: $file (${elapsed}s) → moved ${moved} output files to $DONE_DIR"
    echo "$(date '+%Y-%m-%d %H:%M:%S')  SUCCESS  $file  (${elapsed}s)" >> "$TIMEOUT_LOG"
done

echo "All done."
