# GBWTGraph

GBWTGraph is a [handle graph](https://github.com/vgteam/libhandlegraph) based on the [GBWT](https://github.com/jltsiren/gbwt). Its data model is based on the graph as an alignment of haplotypes. The `gfa2gbwt` tool can be used for converting between a subset of [GFA1](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) and plain and compressed representations of the GBWTGraph.

See [the wiki](https://github.com/jltsiren/gbwtgraph/wiki) for further documentation.

## Overview

The GBWTGraph represents the graph induced by the haplotypes stored in a GBWT index. It uses the GBWT index for graph topology and stores the node sequences in plain form for fast extraction. Construction extracts the sequences from another graph implementing `handlegraph::HandleGraph` or from `gbwtgraph::SequenceSource`.

GBWTGraph supports `handlegraph::HandleGraph` and `handlegraph::SerializableHandleGraph` interfaces. Compared to other handle graph implementations, sequence access is very fast, while graph navigation may be slower. There are also some additional operations:

* `get_sequence_view()` provides direct access to node sequences without decompression, reverse complementation, or memory allocation.
* `follow_paths()` is an analogue of `follow_edges()` using GBWT search states instead of handles. It only follows edges if the resulting path is supported by the haplotypes in the index.
* `compress()` and `decompress()` offer a more space-efficient serialization alternative.

Accessing and decompressing GBWT node records is somewhat slow. Algorithms that repeatedly access the edges in a small subgraph may create a `CachedGBWT` cache using `get_cache()` and pass it explicitly to the relevant queries. Alternatively, they can create a `CachedGBWTGraph` overlay graph that uses a cache automatically. Both types of caches store all accessed records, so a new cache should be created for each subgraph.

GBWTGraph also supports an experimental `SegmentHandleGraph` interface with GFA-like semantics. Each GFA segment with a string name maps to a range of node ids, and GFA links correspond to edges that connect the ends of segments. This interface is currently only available in graphs built using `SequenceSource`.

The package also includes:

* GBWT / GBWTGraph construction from a subset of GFA1, and GFA extraction from a GBWTGraph.
* A minimizer index implementation for indexing the haplotypes in the GBWTGraph.
* GBWT construction from a greedy maximum path cover:
  * Artificial paths that try to cover all length-k contexts equally, either in the entire graph or only in components that do not already contain paths.
  * Concatenations of local length-k haplotypes sampled according to their true frequencies.

## Construction from GFA

The `gfa2gbwt` tool can be used for building GBWTGraph from GFA1, for extracting GFA from the graph, and for converting between plain and compressed representations of GBWTGraph. The tool interprets the GFA file in the following way:

* Overlaps, containments, and tags are ignored.
* Links are induced by the paths; the tool ignores L-lines.
* Experimental W-lines are the primary representation of haplotype paths.
* If there are both P-lines and W-lines in the file, the P-lines are assumed to be reference paths. They are stored with sample name `_gbwt_ref` and with the path name as contig name.
* If there are only P-lines in the file, GBWT metadata can be parsed by providing a regex and a mapping from submatches to metadata fields.

In the plain representation, the GBWT index and the GBWTGraph are stored in separate `.gbwt` and `.gg` files. The compressed representation uses a single `.gbz` file, with the graph stored more space-efficiently than the in-memory representation.

## Dependencies

* [libhandlegraph](https://github.com/vgteam/libhandlegraph) for the handle graph interface.
* [GBWT](https://github.com/jltsiren/gbwt) (latest master) for the backend.
* [SDSL](https://github.com/vgteam/sdsl-lite) (vgteam fork) for low-level data structures.

These dependencies should be installed separately. Because libhandlegraph and SDSL are header-based libraries, having multiple versions of them in the same project may cause issues. Hence all submodules of the main project should use the same copies of these libraries.

All dependencies should be installed before compiling GBWTGraph. By default, libhandlegraph installs to the system directories, while GBWT and SDSL install to the user's home directory.

## Compiling GBWTGraph

GBWTGraph uses C++14 and OpenMP. At the moment, it compiles with g++ (version 6.1 or newer should be enough) on both Mac and Linux. Apple Clang should also work on Mac, but you must install libomp separately from Macports or Homebrew.

Like GBWT, GBWTGraph takes its compiler options from SDSL. For this purpose, you must set `SDSL_DIR` in the makefile to your SDSL main directory before compiling (the default value is `../sdsl-lite`). After that, `make` will compile the library, while `install.sh` will compile and install the headers and the library to your home directory. Another install directory can be specified as `install.sh prefix`.

## CMake build

The CMake build option is provided as a best effort to support some external projects. Instead of using separately installed dependencies, this approach clones GBWT and libhandlegraph as submodules and uses the SDSL from GBWT. As I do not use CMake myself, the build may not always work correctly. If you depend on it, be prepared to fix any issues.
