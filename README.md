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

Run it from the repository root, passing the input CSV and, optionally, an output path:

```bash
./build/main data/input/tree.csv
./build/main data/input/tree.csv path/to/out.csv
```

By default the result is written to `data/output/results.csv`. The output directory is created if it doesn't exist. The default path is relative to the working directory, so run from the repository root. `--help` prints the usage.

## Test

Unit tests for the CSV reader and writer use GoogleTest. It is used from the system if installed, otherwise downloaded by CMake on first configure (needs network). They are built by default:

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
