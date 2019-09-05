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

There is also a minimizer index implementation for indexing the haplotypes in the GBWTGraph.

## Dependencies

* [libhandlegraph](https://github.com/vgteam/libhandlegraph) for the handle graph interface.
* [GBWT](https://github.com/jltsiren/gbwt) (1.0 or later) for the backend.
* [SDSL](https://github.com/simongog/sdsl-lite) for low-level data structures.

These dependencies should be installed separately. Because libhandlegraph and SDSL are header-based libraries, having multiple versions of them in the same project may cause issues. Hence all submodules of the main project should use the same copies of these libraries.

All dependencies should be installed before compiling GBWTGraph. By default, libhandlegraph installs to the system directories, while GBWT and SDSL install to the user's home directory.

## Compiling GBWTGraph

GBWTGraph uses C++11 and OpenMP. At the moment, it compiles with g++ (version 4.9 or newer should be enough) on both Mac and Linux. Clang 9.1 or newer should also work on Mac, but you must install libomp separately from Macports or Homebrew.

Like GBWT, GBWTGraph takes its compiler options from SDSL. For this purpose, you must set `SDSL_DIR` in the makefile to your SDSL main directory before compiling (the default value is `../sdsl-lite`). After that, `make` will compile the library, while `install.sh` will compile and install the headers and the library to your home directory. Another install directory can be specified as `install.sh prefix`.

GBWTGraph is compiled with `-DNDEBUG` by default. Using this option is highly recommended. There are several cases, where SDSL code works correctly but the assertions are incorrect. As SDSL 2.0 is no longer actively supported, we have to wait until the release of SDSL 3.0 to fix these issues.
