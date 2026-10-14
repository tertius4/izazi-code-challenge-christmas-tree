# Christmas Tree Solution

A C++20 command-line tool for the Kaggle Santa 2025 Christmas tree packing challenge. It reads a CSV of tree vertices and, for every N = 1..200, packs N trees (free rotation, no overlap) into the smallest square it can find. The result is a Kaggle submission CSV.

## Setup

Requirements:

- A C++20 compiler (GCC 10+ or Clang 12+)
- CMake 3.16+

On Debian/Ubuntu:

```bash
sudo apt install build-essential cmake
```

## Install / Build

From the repository root:

```bash
cmake -S . -B build
cmake --build build
```

This produces the `build/main` executable. The default build type is Release.

## Use

Run it from the repository root, passing the input CSV and, optionally, an output path:

```bash
./build/main data/input/tree.csv
./build/main data/input/tree.csv path/to/out.csv
```

By default the result is written to `data/output/results.csv` in Kaggle format (`id,x,y,deg`, values prefixed with `s`). A `report.csv` next to it lists, for each N, the square side, its share of the score, and which strategy produced it. The program prints the total score (sum of side²/N; lower is better).

A full run uses all CPU cores (OpenMP) and takes a few minutes. Effort and the strategy cutoffs are set in `SolverConfig` (`src/packing/solver.h`, adjusted in `src/main.cpp`). The output directory is created if it doesn't exist. The default path is relative to the working directory, so run from the repository root. `--help` prints the usage.

## How it works

- **Geometry** (`src/geometry`): rotate and shift a polygon; collision test via edge intersection + point-in-polygon. Touching counts as a collision, so accepted layouts always have a tiny gap.
- **Starting layouts** (`src/packing/starts.*`, `lattice.*`): an exact rotation scan for N=1, random non-overlapping layouts, and interlocking lattices (alternating upright/upside-down columns) at several angles.
- **Simulated annealing** (`src/packing/annealer.*`): moves, rotates or flips one tree at a time, rejecting any move that causes a collision, to shrink the square.
- **Strategy bands** (`src/packing/solver.*`): small N use many random starts, medium N use all lattices plus random starts, large N use the best lattices. Afterwards, a layout for N+1 with one tree removed replaces N if it is smaller.

## Check a result

`scripts/check_submission.py` re-checks a submission independently with shapely (the library Kaggle uses): tree counts, overlaps, coordinate range and the score.

```bash
python3 scripts/check_submission.py data/input/tree.csv data/output/results.csv
```

It needs shapely (`sudo apt install python3-shapely`, or `pip install shapely`).

## Test

Unit tests use GoogleTest. It is used from the system if installed, otherwise downloaded by CMake on first configure (needs network). They are built by default:

```bash
ctest --test-dir build --output-on-failure
```

Pass `-DPACKER_BUILD_TESTS=OFF` to `cmake` to skip them.

## Input format

The input CSV has the columns `vertex,x,y,description`:

```csv
vertex,x,y,description
1,0.0,0.8,Tip
2,0.125,0.5,Right top tier outer
```
