# Christmas Tree Solution

A C++20 command-line tool.
The total score is ±75. Design notes and thought process are in [APPROACH.md](APPROACH.md).

## Quick start

Requirements: a C++20 compiler (GCC 10+ / Clang 12+) and CMake 3.16+ (`sudo apt install build-essential cmake`). OpenMP is used if found. The evaluator also needs `python3` with `pandas` and `shapely`.

```bash
./run.sh            # no args: interactive menu
./run.sh build      # configure + build (Release) -> build/main
./run.sh run        # pack data/input/tree.csv -> data/output/results.csv
./run.sh test       # build + unit tests
./run.sh evaluate   # validate and score the submission with shapely
./run.sh visualise 10   # plot the 10-tree layout
./run.sh visualise-all  # zoomable HTML grid of all 200 layouts (data/output/grid.html)
```

A full run uses all cores and takes about 2-5 minutes. Run from the repository root (the default output path is relative).

## Output

- `data/output/results.csv`: `id,x,y,deg`, values prefixed with `s`.
- `report.csv` next to it: per N, the square side, its share of the score and the strategy that produced it.

The total score is printed at the end. Effort and strategy cutoffs live in `SolverConfig` ([src/packing/solver.h](src/packing/solver.h)), set in [src/main.cpp](src/main.cpp).

## How it works

- **Geometry**: rotate/shift a polygon; collisions via edge intersection + point-in-polygon. Touching counts as a collision, so layouts keep a tiny gap.
- **Starting layouts**: exact rotation scan for N=1, random non-overlapping layouts, and cuts from double lattices (two interlocking trees per cell, the second turned 180°; tilt, spacing and offsets are tuned to make the cell small).
- **Simulated annealing**: moves, rotates or flips one tree at a time, rejecting collisions, plus compaction moves on the trees that set the square's side.
- **Strategy bands**: small N use many random starts, medium N lattices plus random starts, large N the best two lattices. Finally, a layout for N+1 with one tree removed replaces N if smaller.

## Layout

| Path | Contents |
| --- | --- |
| `src/cli`, `src/io` | argument parsing; CSV reading/writing, submission and report output |
| `src/geometry` | polygon and collision code |
| `src/packing` | layout, starts, lattice, annealer, solver |
| `tests/` | GoogleTest unit tests |
| `scripts/evaluator.py` | independent validation, scoring and visualisation |
| `scripts/visualise_grid.py` | zoomable HTML grid of all 200 layouts |
| `data/` | input tree shape and generated output |

## Input format

`vertex,x,y,description`. The shape must be symmetric about the y-axis.

```csv
vertex,x,y,description
1,0.0,0.8,Tip
2,0.125,0.5,Right top tier outer
```

## Tests

```bash
./run.sh test
```

GoogleTest is taken from the system if installed, otherwise downloaded by CMake on first configure. Pass `-DPACKER_BUILD_TESTS=OFF` to skip.
