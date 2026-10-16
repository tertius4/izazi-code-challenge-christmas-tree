# Christmas Tree Solution

A C++20 command-line tool that packs 1 to 200 Christmas trees into the smallest squares.
Runs are incremental: each one starts from the best layouts found so far and only saves
improvements, so the score keeps getting better the more you run it. The current score is the
`total` row of `data/current-best/report.csv`. Design notes and thought process are in
[APPROACH.md](APPROACH.md).

## Quick start

Requirements: a C++20 compiler (GCC 10+ / Clang 12+) and CMake 3.16+ (`sudo apt install build-essential cmake`). OpenMP is used if found. The evaluator also needs `python3` with `pandas` and `shapely`.

```bash
./run.sh            # no args: interactive menu
./run.sh build      # configure + build (Release) -> build/main
./run.sh run        # improve data/current-best for 10 minutes on half the cores
./run.sh run --minutes 60 --threads 4   # longer, lighter on the machine
./run.sh run --n 13-30                  # only work on some n
./run.sh test       # build + unit tests
./run.sh evaluate   # validate and score data/current-best/results.csv with shapely
./run.sh visualise 10   # plot the 10-tree layout
./run.sh visualise-all  # zoomable HTML grid of all 200 layouts (grid.html next to results.csv)
```

Run from the repository root (the default paths are relative). Every improvement is printed as
it happens. Ctrl-C stops a run early and still saves. The run saves every two minutes anyway,
writing to a temporary file first, so even a killed run never damages the store.

If `data/current-best/results.csv` doesn't exist, the first run builds every layout from scratch
(lattice cuts, random starts, neighbour passes; this takes a long time) and then improves it.

## Output

`data/current-best/` (pass `--best DIR` to use another folder):
- `results.csv`: `id,x,y,deg`, values prefixed with `s`.
- `report.csv`: per N, the square side, its share of the score and how the layout was last
  improved.

## How it works

- **The tree**: hard-coded, not read from a file. It is exactly the union of four convex parts
  (the three tiers, each extended to a triangle, and the trunk). Two trees collide when any
  pair of parts does, which a separating-axis test on 3–4 vertices decides quickly. Touching
  counts as a collision, so layouts keep a tiny gap. Tests check this against the general
  polygon test on a million random pairs and on pairs right at the point of contact.
- **Layouts** keep a grid over the trees, so a collision check only looks at nearby trees.
- **Simulated annealing**: moves, rotates or flips one tree at a time, rejecting collisions,
  plus compaction moves on the trees that set the square's side.
- **Improving** (`src/packing/improver.cpp`): passes go through n from the largest down to
  the smallest, threads sharing the pass. Each visit of n tries, in order: n+1 trees minus each
  of its best few trees to remove; n-1 trees plus one, squeezed; squeezing the best layout
  until it fails twice in a row; a short re-anneal; and for small n, every other pass, a fresh
  random start. A result replaces the best only if it is valid and smaller. An improvement
  makes n+1 try "n plus a tree" right away.
- **Squeezing** (`src/packing/squeezer.cpp`, by far the most productive): pull the trees into a
  square slightly smaller than the best, letting them overlap, then anneal single-tree moves
  that shrink the overlap (the penetration depth between the trees' convex parts) until none is
  left. How much to shrink adapts per n: doubled after a success, halved after a failure.
  Allowing overlap gives the search a slope to follow; the plain annealer, which rejects every
  colliding move, finds no move that shrinks a tight packing.
- **Building from scratch** (`src/packing/solver.cpp`, only without a store): an exact rotation
  scan for N=1, random starts for small N, and cuts from double lattices (two interlocking trees
  per cell, the second turned 180°) for larger N.

## Layout

| Path | Contents |
| --- | --- |
| `src/cli`, `src/io` | argument parsing; CSV reading/writing, submission and report output |
| `src/tree` | the tree: outline, convex parts, collision test |
| `src/geometry` | general polygon code (used by tests, lattice cuts and rotation) |
| `src/packing` | layout, starts, lattice, annealer, rotation, improver, solver |
| `tests/` | GoogleTest unit tests |
| `scripts/evaluator.py` | independent validation, scoring and visualisation |
| `scripts/visualise_grid.py` | zoomable HTML grid of all 200 layouts |
| `data/` | the tree outline (used by tests) and the best-known layouts |

## Tests

```bash
./run.sh test
```

GoogleTest is taken from the system if installed, otherwise downloaded by CMake on first configure. Pass `-DPACKER_BUILD_TESTS=OFF` to skip.
