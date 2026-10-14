# Christmas Tree Solution

A C++20 command-line tool that reads a CSV of tree vertices and writes a CSV result. The packing logic is not implemented yet, so for now the input is copied to the output.

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

Run it from the repository root, passing the input CSV:

```bash
./build/main data/input/tree.csv
```

The result is written to `data/output/result.csv`. The output directory is created if it doesn't exist. Run from the repository root, because the output path is relative to the working directory.

The input CSV has the columns `vertex,x,y,description`:

```csv
vertex,x,y,description
1,0.0,0.8,Tip
2,0.125,0.5,Right top tier outer
```
