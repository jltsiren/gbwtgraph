# GBWTGraph

This is a [handle graph](https://github.com/vgteam/libhandlegraph) implementation based on the [GBWT](https://github.com/jltsiren/gbwt). The development started as a part of [VG](https://github.com/vgteam/vg), but the implementation was moved to an independent library with minimal dependencies.

See [the wiki](https://github.com/jltsiren/gbwtgraph/wiki) for further documentation.

## Overview

The GBWTGraph represents the graph induced by the haplotypes stored in a GBWT index. It uses the GBWT index for graph topology and stores the node sequences in plain form for fast extraction. The construction requires a sequence source, which must implement the following subset of `handlegraph::HandleGraph` interface for all nodes in forward orientation:

* `get_handle()`
* `get_length()`
* `get_sequence()`

GBWTGraph supports `handlegraph::HandleGraph` and `handlegraph::SerializableHandleGraph` interfaces. Compared to other handle graph implementations, sequence access is very fast, while graph navigation may be slower. There are also some additional operations:

* `get_sequence_view()` provides direct access to node sequences without decompression, reverse complementation, or memory allocation.
* `follow_paths()` is an analogue of `follow_edges()` using GBWT search states instead of handles. It only follows edges if the resulting path is supported by the haplotypes in the index.

Accessing and decompressing GBWT node records is somewhat slow. Algorithms that repeatedly access the edges in a small subgraph may create a `CachedGBWT` cache using `get_cache()` and pass it explicitly to the relevant queries. Alternatively, they can create a `CachedGBWTGraph` overlay graph that uses a cache automatically. Both types of caches store all accessed records, so a new cache should be created for each subgraph.

The package also includes:

* Direct GBWT / GBWTGraph construction from GFA1 (assuming integer segment identifiers, no overlaps, and no containments).
* A minimizer index implementation for indexing the haplotypes in the GBWTGraph.
* GBWT construction from a greedy maximum path cover:
  * Artificial paths that try to cover all length-k contexts equally, either in the entire graph or only in components that do not already contain paths.
  * Concatenations of local length-k haplotypes sampled according to their true frequencies.

## Dependencies

* [libhandlegraph](https://github.com/vgteam/libhandlegraph) for the handle graph interface.
* [GBWT](https://github.com/jltsiren/gbwt) (latest master) for the backend.
* [SDSL](https://github.com/vgteam/sdsl-lite) (vgteam fork) for low-level data structures.

These dependencies should be installed separately. Because libhandlegraph and SDSL are header-based libraries, having multiple versions of them in the same project may cause issues. Hence all submodules of the main project should use the same copies of these libraries.

All dependencies should be installed before compiling GBWTGraph. By default, libhandlegraph installs to the system directories, while GBWT and SDSL install to the user's home directory.

## Compiling GBWTGraph

GBWTGraph uses C++14 and OpenMP. At the moment, it compiles with g++ (version 6.1 or newer should be enough) on both Mac and Linux. Apple Clang  should also work on Mac, but you must install libomp separately from Macports or Homebrew.

Like GBWT, GBWTGraph takes its compiler options from SDSL. For this purpose, you must set `SDSL_DIR` in the makefile to your SDSL main directory before compiling (the default value is `../sdsl-lite`). After that, `make` will compile the library, while `install.sh` will compile and install the headers and the library to your home directory. Another install directory can be specified as `install.sh prefix`.

## CMake build

There is an alternative build option using CMake. Instead of using separately installed dependencies, this approach clones them as submodules.
