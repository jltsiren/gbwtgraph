# Simple-SDS serialization format

Version 3. Updated 2021-05-12.

## Basics

This document specifies the portable simple-sds serialization format for GBWTGraph.
It is an alternative to the old file format based on SDSL data structures.
The format builds upon the data structures described in:

* <https://github.com/jltsiren/simple-sds/blob/main/SERIALIZATION.md>
* <https://github.com/jltsiren/gbwt/blob/master/SERIALIZATION.md>

## GBWTGraph

**GBWTGraph** represents a bidirected sequence graph induced by a set of paths.
It uses a GBWT index for graph topology and for storing the paths.
The GBWT must be bidirectional and contain metadata with path names.
Paths with sample name `_gbwt_ref` are assumed to be reference paths.

Serialization format for GBWTGraph:

1. GBWTGraph header
2. Tags
3. Sequences
4. Node-to-segment translation

**Tags** are stored in the same format as in the GBWT.
Key `source` indicates the GBWTGraph implementation used for serializing the graph.
The original implementation corresponds to value `jltsiren/gbwtgraph`.

### GBWTGraph header

**GBWTGraph header** is a 24-byte (3-element) structure with the following fields:

1. `tag`: Value `0x6B3764AF` as a 32-bit integer.
2. `version`: File format version as a 32-bit integer.
3. `nodes`: Number of nodes in the graph as an element.
4. `flags`: Binary flags as an element.

The first two fields are 32-bit unsigned little-endian integers for compatibility with the SDSL-based serialization format.
Simple-SDS serialization format requires file format version `3`.

Field `nodes` counts the number of original nodes **present** in the graph.
A node is present if the local alphabet size in the corresponding GBWT nodes is nonzero.
If the GBWT stores a subgraph induced by the paths, this is equivalent to requiring that the original node is visited by a path.

The following flags are supported:

* `0x0001`: The translation structure is present.
* `0x0002`: Simple-SDS format.

Other flag bits must not be set.
The translation flag must be set if and only if the translation structure is nonempty.
If the simple-sds bit is not set, the serialized data is in the SDSL format.

### Sequences

The **sequences** structure attaches a string label to every original node present in the graph.

Serialization format for sequences:

1. `sequences`: Node labels as a string array.

The label of original node `v` is string `v - floor(offset / 2) - 1` in the string array, where `offset` is the alphabet offset in the GBWT.
This usually means that the nodes map to consecutive identifiers starting from `0`.
If a node is not present in the graph, the corresponding string must be empty.

### Node-to-segment translation

The GBWT uses nodes with integer identifiers, while the [GFA specification](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) uses **segments** with string names.
GBWTGraph needs to store a **translation** between node identifiers and string names if the names cannot be interpreted as integers.
Because some applications expect nodes with short sequences, the translation also supports mapping segments to ranges of nodes.

Serialization format for translation:

1. `segments`: Segment names as a string array.
2. `mapping`: First node identifier mapping to each segment as a sparse vector.

If the translation is not in use, both structures must be empty.
Otherwise all original nodes in the closed interval from `1` to `nodes` must be present.
The segment with the `i`th string as its name then corresponds to the concatenation of nodes from `mapping.select(i)` (inclusive) to `mapping.select(i + 1)` (exclusive).

**Note:** Because the translation is a core part of the GBWTGraph, it is not serialized as an optional structure.
