# GBZ file format

GBZ version 1, GBWTGraph version 3. Updated 2021-10-11.

## Basics

This document specifies the GBZ file format for storing GFA with many paths.
It includes a portable simple-sds serialization format for GBWTGraph as an alternative to the old SDSL-based format.
The format builds upon the data structures described in:

* <https://github.com/jltsiren/simple-sds/blob/main/SERIALIZATION.md>
* <https://github.com/jltsiren/gbwt/blob/master/SERIALIZATION.md>

The notation uses simple-sds conventions.
Array indexes start at 0.
Rank queries count the number of occurrences strictly before the query position, as in SDSL.
A rank query past the end of the sequence counts the total number of occurrences in the sequence.
Select queries use 0-based indexing (SDSL uses 1-based indexing).
A past-the-end select query should be understood as the length of the sequence, though this is usually clarified later.

## GBZ

**GBZ** is a space-efficient binary format for a subset of [GFA1](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) with many paths.
It ignores optional GFA fields and assumes that there are no overlaps or containments.
The format builds upon the simple-sds formats for serializing the GBWT and the GBWTGraph.

GBZ file format:

1. GBZ header
2. Tags
3. GBWT
4. GBWTGraph

**Tags** are stored in the same format as in the GBWT.
Key `source` indicates the implementation that generated the file.
The original implementation corresponds to value `jltsiren/gbwtgraph`.

### GBZ header

**GBZ header** is a 16-byte (2-element) structure with the following fields:

1. `tag`: Value `0x205A4247` as a 32-bit integer (corresponds to string `GBZ `).
2. `version`: File format version as a 32-bit integer.
3. `flags`: Binary flags as an element.

The first two fields are 32-bit unsigned little-endian integers for consistency with the other headers.
Current file format version is 1.
Flags are not used in version 1.

### GBWT / GFA paths

If the GBWT does not contain metadata with path names, each original path in the GBWT corresponds to a P-line in GFA.
The integer identifier of each original path is used as the `PathName` field.

Otherwise the original paths for sample `_gbwt_ref` are assumed to be reference paths that correspond to GFA P-lines.
Each reference path must have a different contig identifier in the GBWT path name.
The corresponding contig name is used as the `PathName` field.

Paths for all other samples are assumed to be haplotypes that correspond to GFA W-lines.
GBWT path name fields map to GFA W-line fields in the following way:

* The sample name corresponding to `sample` maps to `SampleId`.
* The contig name corresponding to `contig` maps to `SeqId`.
* `phase` maps to `HapIndex`.
* `fragment` maps to `SeqStart` (`SeqEnd` can be derived from `fragment` and the path itself).

If the GBWT metadata does not contain sample (contig) names, the integer identifier of the sample (contig) is used instead.

## GBWTGraph

**GBWTGraph** represents a bidirected sequence graph induced by a set of paths.
It uses a GBWT index for graph topology and for storing the paths.
The GBWT must be bidirectional.

Serialization format for GBWTGraph:

1. GBWTGraph header
2. Sequences
3. Node-to-segment translation

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
As the GBWT stores a subgraph induced by the paths, this is equivalent to requiring that the original node is visited by a path.

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

**Note:** If a node is not present in the graph, the corresponding string should be empty.
While this is not a strict requirement, it is a good practice that makes the serialization deterministic.

### Node-to-segment translation

The GBWT uses nodes with integer identifiers, while the [GFA specification](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) uses **segments** with string names.
GBWTGraph needs to store a **translation** between node identifiers and string names if the names cannot be interpreted as integers.
Because some applications expect nodes with short sequences, the translation also supports mapping segments to ranges of nodes.

Serialization format for translation:

1. `segments`: Segment names as a string array.
2. `mapping`: First node identifier mapping to each segment as a sparse vector.

If the translation is not in use, both structures must be empty.
Otherwise the original nodes must have consecutive identifiers starting from `1`.
The segment with the `i`th string as its name then corresponds to the concatenation of nodes from `mapping.select(i)` (inclusive) to `mapping.select(i + 1)` (exclusive).

**Note:** Because the translation is a core part of the GBWTGraph, it is not serialized as an optional structure.

**Note:** For the last segment, the end of the node interval is the length of `mapping`.

**Note:** Some segments may consist of nodes that are not present in the graph.
As with sequences, such segments should have empty names to make the serialization deterministic.
