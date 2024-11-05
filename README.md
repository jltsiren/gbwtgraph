# GBWTGraph

GBWTGraph is a [handle graph](https://github.com/vgteam/libhandlegraph) based on the [GBWT](https://github.com/jltsiren/gbwt). Its data model is based on the graph as an alignment of haplotypes. The `gfa2gbwt` tool can be used for converting between a subset of [GFA1](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md), the GBZ file format, and GBWTGraph.

See [the wiki](https://github.com/jltsiren/gbwtgraph/wiki) for further documentation.

There is also a partial [Rust implementation](https://github.com/jltsiren/gbwt-rs) of the GBWTGraph.

## Overview

The GBWTGraph represents the graph induced by the haplotypes stored in a GBWT index. It uses the GBWT index for graph topology and stores the node sequences in plain form for fast extraction. Construction extracts the sequences from another graph implementing `handlegraph::HandleGraph` or from `gbwtgraph::SequenceSource`.

GBWTGraph supports the following libhandlegraph interfaces:

* `HandleGraph`
* `PathHandleGraph`
* `NamedNodeBackTranslation`
* `SerializableHandleGraph` (the GBWT index must be serialized/deserialized separately)

Compared to other handle graph implementations, sequence access is very fast, while graph navigation may be slower. There are also some additional operations:

* `get_sequence_view()` provides direct access to node sequences without decompression, reverse complementation, or memory allocation.
* `follow_paths()` is an analogue of `follow_edges()` using GBWT search states instead of handles. It only follows edges if the resulting path is supported by the haplotypes in the index.
* `simple_sds_serialize()` and `simple_sds_load()` offer a more space-efficient serialization alternative.

Accessing and decompressing GBWT node records is somewhat slow. Algorithms that repeatedly access the edges in a small subgraph may create a `CachedGBWT` cache using `get_cache()` and pass it explicitly to the relevant queries. Alternatively, they can create a `CachedGBWTGraph` overlay graph that uses a cache automatically. Both types of caches store all accessed records, so a new cache should be created for each subgraph.

GBWTGraph also supports an experimental `SegmentHandleGraph` interface with GFA-like semantics. Each GFA segment with a string name maps to a range of node ids, and GFA links correspond to edges that connect the ends of segments. This interface is currently only available in graphs built using `SequenceSource`.

The package also includes:

* [GBZ file format](https://github.com/jltsiren/gbwtgraph/blob/master/SERIALIZATION.md).
* GBWT / GBWTGraph construction from a subset of GFA1, and GFA extraction from a GBWTGraph.
* A generic kmer index and a generic minimizer index for indexing the haplotypes in the GBWTGraph.
* GBWT construction from a greedy maximum path cover:
  * Artificial paths that try to cover all length-k contexts equally, either in the entire graph or only in components that do not already contain paths.
  * Concatenations of local length-k haplotypes sampled according to their true frequencies.
* Support for paralellizing GBWT construction over weakly connected components of the graph:
  * `gbwt_construction_jobs()` for determining the construction jobs.
  * `assign_paths()` and `insert_paths()` for passing reference paths to the new GBWT.
  * `MetadataBuilder` for generating GBWT metadata.

## Construction from GFA

The `gfa2gbwt` tool can be used for building GBWTGraph from GFA1, for extracting GFA from the graph, and for converting between plain and compressed representations of GBWTGraph. The tool interprets the GFA file in the following way:

* Overlaps, containments, and tags are ignored.
* Segments and links are inferred from the paths; the corresponding S-lines and L-lines must still exist in the file.
* Experimental W-lines are the primary representation of haplotype paths.
* If there are both P-lines and W-lines in the file, the P-lines are assumed to be reference paths. They are stored with sample name `_gbwt_ref` and with the path name as contig name.
* If there are only P-lines in the file, GBWT metadata can be parsed by providing a regex and a mapping from submatches to metadata fields.

In the plain representation, the GBWT index and the GBWTGraph are stored in separate `.gbwt` and `.gg` files. The compressed representation uses a single `.gbz` file, with the graph stored more space-efficiently than the in-memory representation.

## Sequence extraction

Sequences corresponding to paths can be extracted from a GBZ graph using the `gbz_extract` tool. The tool extracts sequences using multiple threads and writes them to the standard output. Each sequence is terminated by an endline character (`'\n`), making this tool suitable for multi-string BWT construction using tools such as [grlBWT](https://github.com/ddiazdom/grlBWT).

The extracted sequences are written in the same order the paths appear in the GBWT. By default, the line number in the output (which becomes the sequence identifier in the multi-string BWT) is the same as the path identifier in the GBWT. If option `-b` / `--both-orientations` is used, each sequence is followed by the reverse complement of the same sequence. This can be used for building an FMD-index. When reverse complements are included, the line number is the same as the GBWT sequence identifier.

Instead of extracting all sequences, it is possible to limit the extraction to the weakly connected component(s) corresponding to a contig name (e.g. `chr19`). This can be done using option `-c` / `--contig`. When this option is used, the sequences will maintain the same relative order as the corresponding paths in the GBWT.

The number of extraction threads can be changed using option `-t` / `--threads`. In addition to the extraction threads, there will be a main thread responsible for writing the sequences to the output. The extraction process is often I/O bound, and it is unlikely that more than a few threads will be active at the same time.

## Dependencies

* [libhandlegraph](https://github.com/vgteam/libhandlegraph) for the handle graph interface.
* [GBWT](https://github.com/jltsiren/gbwt) for the backend.
* [SDSL](https://github.com/vgteam/sdsl-lite) (vgteam fork) for low-level data structures.

These dependencies should be installed separately (the latest master should always work). Because libhandlegraph and SDSL are header-based libraries, having multiple versions of them in the same project may cause issues. Hence all submodules of the main project should use the same copies of these libraries.

All dependencies should be installed before compiling GBWTGraph. By default, libhandlegraph installs to system directories, while GBWT and SDSL install to the user's home directory. Dependencies not installed in system directories should use the same install prefix as SDSL.

## Compiling GBWTGraph

This library is designed to take the compiler options from the [vgteam fork](https://github.com/vgteam/sdsl-lite) of the Succinct Data Structures Library 2.0 (SDSL). It currently requires a recent C++ compiler supporting C++17 and OpenMP. GCC is recommended, as the multithreaded `std::sort` from libstdc++ parallel mode speeds up some algorithms. On Apple systems, GBWTGraph can be built with Apple Clang, but libomp must be installed via Macports or Homebrew.

GBWTGraph is frequently tested in the following environments:

* Intel Linux (Ubuntu) with GCC.
* Intel macOS with GCC and Apple Clang.
* ARM macOS with Apple Clang.

Before compiling, you must set `SDSL_DIR` in the makefile to your SDSL main directory. The default value is `../sdsl-lite`, which is usually appropriate. The makefile will read `$SDSL_DIR/Make.helper` to determine compilers and compiler options.

After that, `make` will compile the library, while `install.sh` will compile and install the headers and the library to your home directory. Another install directory can be specified with `install.sh prefix`.

## Citing GBWTGraph

Jouni Sirén, Jean Monlong, Xian Chang, Adam M. Novak, Jordan M. Eizenga, Charles Markello, Jonas A. Sibbesen, Glenn Hickey, Pi-Chuan Chang, Andrew Carroll, Namrata Gupta, Stacey Gabriel, Thomas W. Blackwell, Aakrosh Ratan, Kent D. Taylor, Stephen S. Rich, Jerome I. Rotter, David Haussler, Erik Garrison, and Benedict Paten:
**Pangenomics enables genotyping of known structural variants in 5202 diverse genomes**.
Science 374(6574):abg8871, 2021.
DOI: [10.1126/science.abg8871](https://doi.org/10.1126/science.abg8871)

## Citing GBZ File Format

Jouni Sirén and Benedict Paten:
**GBZ file format for pangenome graphs**.
Bioinformatics 38(22):5012-5018, 2022.
DOI: [10.1093/bioinformatics/btac656](https://doi.org/10.1093/bioinformatics/btac656)
