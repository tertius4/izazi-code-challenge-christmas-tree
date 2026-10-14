#!/usr/bin/env bash
# Simple CLI for building, running and testing the packer.
set -euo pipefail

cd "$(dirname "$0")"

BUILD_DIR=build
INPUT=data/input/tree.csv
OUTPUT=data/output/results.csv

usage() {
    cat <<USAGE
Usage: ./run.sh <command> [args]

Commands:
  build [Release|Debug]   Configure and build (default: Release)
  run [input] [output]    Build, then run the packer (default input: $INPUT)
  test                    Build, then run the unit tests
  evaluate [submission]   Validate and score a submission (default: $OUTPUT)
  visualise N [submission]  Plot the N-tree layout from a submission
  clean                   Remove the build directory
  menu                    Interactive numbered menu (default with no command)
  help                    Show this message
USAGE
}

build() {
    cmake -S . -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE="${1:-Release}"
    cmake --build "$BUILD_DIR" -j"$(nproc)"
}

menu() {
    local items=("Build" "Run" "Test" "Evaluate" "Visualise" "Clean" "Help" "Quit")
    while true; do
        echo
        for i in "${!items[@]}"; do
            echo "  $((i + 1))) ${items[$i]}"
        done
        read -rp "Choose [1-${#items[@]}]: " choice || exit 0
        # A failing action shouldn't close the menu.
        set +e
        case "$choice" in
            1) build ;;
            2) read -rp "Input [$INPUT]: " in; read -rp "Output [$OUTPUT]: " out
               "$0" run "${in:-$INPUT}" "${out:-$OUTPUT}" ;;
            3) "$0" test ;;
            4) read -rp "Submission [$OUTPUT]: " sub; "$0" evaluate "${sub:-$OUTPUT}" ;;
            5) read -rp "Tree count N: " n; read -rp "Submission [$OUTPUT]: " sub
               "$0" visualise "$n" "${sub:-$OUTPUT}" ;;
            6) "$0" clean ;;
            7) usage ;;
            8) exit 0 ;;
            *) echo "Invalid choice: $choice" ;;
        esac
        set -e
    done
}

cmd="${1:-menu}"
shift || true

case "$cmd" in
    build) build "$@" ;;
    run)
        in="${1:-$INPUT}"
        shift || true
        ./"$BUILD_DIR"/main "$in" "$@"
        ;;
    test)
        build
        ctest --test-dir "$BUILD_DIR" --output-on-failure
        ;;
    evaluate) python3 scripts/evaluator.py "${1:-$OUTPUT}" ;;
    visualise|visualize)
        [[ $# -ge 1 ]] || { echo "visualise needs a tree count N" >&2; exit 1; }
        n="$1"
        python3 scripts/evaluator.py "${2:-$OUTPUT}" --visualize "$n"
        ;;
    clean) rm -rf "$BUILD_DIR" ;;
    menu) menu ;;
    help|-h|--help) usage ;;
    *) echo "Unknown command: $cmd" >&2; usage >&2; exit 1 ;;
esac
