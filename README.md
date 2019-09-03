# GBWTGraph

This is a [Handle Graph](https://github.com/vgteam/libhandlegraph) implementation based on the [GBWT](https://github.com/jltsiren/gbwt). The development started as a part of [VG](https://github.com/vgteam/vg), but the implementation was moved into an independent library with minimal dependencies.

See the wiki for further documentation.

## Overview

GBWTGraph represents the graph induced by the haplotypes stored in a GBWT index. It uses the GBWT index for graph topology and stores the node sequences in plain form for fast extraction. The construction requires a sequence source, which must implement the following subset of `handlegraph::HandleGraph` interface for all nodes in forward orientation:

* `get_handle()`
* `get_length()`
* `get_sequence()`

In addition to the standard `handlegraph::HandleGraph` and `handlegraph::SerializableHandleGraph` interfaces, GBWTGraph supports a number of additional operations:

* `get_sequence_view()` provides direct access to node sequences without decompression, reverse complementation, or memory allocation.
* `follow_paths()` is an analogue of `follow_edges()` using GBWT search states instead of handles. It only follows edges if the resulting path is supported by the haplotypes in the index.

There is also a minimizer index implementation for indexing the haplotypes in the GBWTGraph.

## Dependencies

* [libhandlegraph](https://github.com/vgteam/libhandlegraph) for the Handle Graph interface.
* [GBWT](https://github.com/jltsiren/gbwt) 1.0 or later for the backend.
* [SDSL](https://github.com/simongog/sdsl-lite) for low-level data structures.

These dependencies should be installed separately. Because libhandlegraph and SDSL are header-based libraries, having multiple versions of them in the same project may cause issues.

## Compiling GBWTGraph
