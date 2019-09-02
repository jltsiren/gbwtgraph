# GBWTGraph

This is a [Handle Graph](https://github.com/vgteam/libhandlegraph) implementation based on the [GBWT](https://github.com/jltsiren/gbwt). The development started as a part of [VG](https://github.com/vgteam/vg), but the implementation was moved into an independent library with minimal dependencies.

See the wiki for further documentation.

## Overview

## Dependencies

* [libhandlegraph](https://github.com/vgteam/libhandlegraph) for the Handle Graph interface.
* [GBWT](https://github.com/jltsiren/gbwt) 1.0 or later for the backend.
* [SDSL](https://github.com/simongog/sdsl-lite) for low-level data structures.

These dependencies should be installed separately. Because libhandlegraph and SDSL are header-based libraries, having multiple versions of them in the same project may cause issues.

## Compiling GBWTGraph
