# Christmas Tree Solution

A C++20 command-line tool that packs 1 to 200 Christmas trees into the smallest squares.
Runs are incremental: each one starts from the best layouts found so far and only saves
improvements, so the score keeps getting better the more you run it. The current score is the
`total` row of `data/current-best/report.csv`. Design notes and thought process are in
[APPROACH.md](APPROACH.md).

## Quick start

Requirements: a C++20 compiler (GCC 10+ / Clang 12+) and CMake 3.16+ (`sudo apt install build-essential cmake`). The evaluator also needs `python3` with `pandas` and `shapely`.

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
- `crystals.csv`: the densest periodic patterns found (see below), one row per tree of each
  pattern's cell. Searched for once, on the first run without it (a few minutes), then reused.

## How it works

- **The tree**: hard-coded, not read from a file. It is exactly the union of four convex parts
  (the three tiers, each extended to a triangle, and the trunk). Two trees collide when any
  pair of parts does, which a separating-axis test on 3–4 vertices decides quickly. Touching
  counts as a collision, so layouts keep a tiny gap. Tests check this against the general
  polygon test on a million random pairs and on pairs right at the point of contact.
- **Layouts** keep a grid over the trees, so a collision check only looks at nearby trees.
- **Simulated annealing**: moves, rotates or flips one tree at a time, rejecting collisions,
  plus compaction moves on the trees that set the square's side.
- **Improving** (`src/packing/improver.cpp`): threads take one n at a time. Half the visits
  go to the n with the best score gain per second so far this run, the rest to a sweep from
  the largest n down. Each visit of n tries, in order: n+1 trees minus each of its best few
  trees to remove; n-1 trees plus one (at the best spot or a random one); squeezing the best
  layout until it fails often enough; a border squeeze (only trees near the walls move);
  compressing until it fails often enough; shaking; crossover; a crystal cut (n ≥ 20); dense
  starts (n ≤ 40, more for n ≤ 10); a walk of re-anneals at rising temperatures, each squeezed
  down; and for small n, every other visit, a fresh random start. A result replaces the best
  only if it is valid and smaller. A real improvement makes n+1 try "n plus a tree" right
  away. Each operator (but the squeeze) runs with a chance set by its score gain per second
  in that band of n, so operators that don't pay off there are mostly skipped. The stats per
  operator and band are printed on every save and written to `stats.txt`.
- **Squeezing** (`src/packing/squeezer.cpp`, by far the most productive): pull the trees into a
  square slightly smaller than the best, letting them overlap, then move every tree at once
  down the slope of the overlap (L-BFGS on the penetration depth between the trees' convex
  parts), with single-tree annealing moves as a fallback. How much to shrink adapts per n:
  doubled after a success, halved after a failure. Allowing overlap gives the search a slope
  to follow; the plain annealer, which rejects every colliding move, finds no move that shrinks
  a tight packing. Refinements:
  - *guided local search*: when the slope runs out with overlap left, the pairs that still
    overlap get a heavier penalty (the deepest most), which pushes the search out of the jam;
  - *line moves*: some single-tree moves try several distances along x, y or the turn and take
    the best;
  - *border only*: the inside of the layout can stay fixed while the trees near the walls move;
  - *symmetric*: a layout built point-symmetric about the square's centre can be kept so,
    which halves the variables.
- **Compressing**: the square's side becomes a variable too: minimise side + weight·overlap²
  with the weight rising, so the walls close in as far as the trees let them; then squeeze
  into that side to make it valid.
- **Shaking**: swap, flip or move a few trees, settle the result back into the best square,
  and walk on from there; each settled layout gets one squeeze into a smaller square.
- **Crossover** (`src/packing/crossover.cpp`): each n keeps a few good, different layouts.
  A child takes the best layout's trees on one side of a random line and another layout's (an
  elite, the best turned or mirrored, or n±1's best) on the other, then drops or adds trees
  to get n and is settled.
- **Adding a tree at the best spot** (`src/packing/nfp.cpp`): the tree's no-fit polygon against
  another (where the two overlap) is the union of one convex polygon per pair of convex parts.
  Its corners are where a new tree nestles against an old one; those and a grid are tried at
  several angles, and the spot with the least overlap wins.
- **Crystals** (`src/packing/crystal.cpp`): periodic patterns with a cell of any shape holding
  2–4 trees, made as dense as possible by annealing the cell. The best is about 0.3043 per
  tree (81% of the area). For each n, the pattern is turned to many angles, cut into the
  smallest square holding n whole trees, compressed (border first) and squeezed down.
- **Dense starts**: random overlapping trees (half the time point-symmetric) in a roomy
  square, compressed, squeezed valid, and squeezed down if they came close to the best.
- **Building from scratch** (`src/packing/solver.cpp`, only without a store): an exact rotation
  scan for N=1, random starts for small N, and cuts from double lattices (two interlocking trees
  per cell, the second turned 180°) for larger N.

## Layout

| Path | Contents |
| --- | --- |
| `src/cli`, `src/io` | argument parsing; CSV reading/writing, submission and report output |
| `src/tree` | the tree: outline, convex parts, collision test |
| `src/geometry` | general polygon code (used by tests, lattice cuts and rotation) |
| `src/util` | `parallelFor`: runs a loop on a pool of `std::thread`s |
| `src/packing` | layout, starts, lattice, crystal, annealer, squeezer, nfp, crossover, rotation, improver, solver |
| `tests/` | GoogleTest unit tests |
| `scripts/evaluator.py` | independent validation, scoring and visualisation |
| `scripts/visualise_grid.py` | zoomable HTML grid of all 200 layouts |
| `data/` | the tree outline (used by tests) and the best-known layouts |

## Tests

```bash
./run.sh test
```

GoogleTest is taken from the system if installed, otherwise downloaded by CMake on first configure. Pass `-DPACKER_BUILD_TESTS=OFF` to skip.
