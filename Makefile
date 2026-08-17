SDSL_DIR=../sdsl-lite
include $(SDSL_DIR)/Make.helper

BUILD_BIN=bin
BUILD_LIB=lib
BUILD_OBJ=obj
SOURCE_DIR=src

# --- Optional build-time behaviors, overridable on the command line (e.g.
# `make GBWTGRAPH_USE_OPENMP=0`). Defaults match GBWTGraph's traditional
# behavior. See include/gbwtgraph/gbwtgraph.h / gbz.h (GBZ, GBWTGraph) for
# what each one changes.

# Build with OpenMP parallelism (1) or without it (0).
GBWTGRAPH_USE_OPENMP=1
# GBZ/GBWTGraph's shared-memory-backed node sequence storage reuses gbwt's
# StringArray, so it only exists when gbwt itself was built with shared
# memory. include/gbwtgraph/utils.h defines the GBWTGRAPH_ENABLE_SHARED_MEMORY
# macro directly from gbwt/config.h (see GBWT_ENABLE_SHARED_MEMORY there), so
# this Makefile is not responsible for passing it and cannot disagree with
# gbwt about it; it only checks the same header, below, to decide whether to
# link librt.
GBWTGRAPH_ENABLE_SHARED_MEMORY=$(shell echo '\#include <gbwt/config.h>' | $(MY_CXX) -I$(INC_DIR) -E -x c++ - 2>/dev/null | grep -q GBWT_ENABLE_SHARED_MEMORY && echo 1 || echo 0)

# gbwt's own headers make the same behavioral choices under matching GBWT_*
# macros, so passing mismatched ones here can produce confusing compile or
# link errors against however ../gbwt was actually built.
#
# For OpenMP, this Makefile passes gbwt's matching macro whenever the
# corresponding GBWTGRAPH_* one is set. Shared memory does not follow this
# pattern; see GBWTGRAPH_ENABLE_SHARED_MEMORY above.
DEFINES=

# Multithreading with OpenMP.
ifeq ($(GBWTGRAPH_USE_OPENMP), 1)
    PARALLEL_FLAGS=-fopenmp -pthread -DGBWTGRAPH_USE_OPENMP -DGBWT_USE_OPENMP
else
    # Without OpenMP, the "#pragma omp ..." directives sprinkled through the
    # source are inert (silently ignored by the compiler), but -Wall turns on
    # -Wunknown-pragmas, which would otherwise warn about every one of them.
    PARALLEL_FLAGS=-pthread -Wno-unknown-pragmas
endif

# Directories for dependencies.
INCLUDES=-Iinclude -I$(INC_DIR)
LIBS=-L$(LIB_DIR) -lgbwt -lhandlegraph -lsdsl

# Use pkg-config to find system dependencies.
ifeq ($(shell pkg-config --exists libcrypto libzstd && echo 1), 1)
    $(info Found libcrypto and libzstd.)
    INCLUDES += $(shell pkg-config --cflags libcrypto libzstd)
    LIBS += $(shell pkg-config --libs libcrypto libzstd)
else
    $(error Could not find libcrypto or libzstd. Please update PKG_CONFIG_PATH)
endif

ifeq ($(GBWTGRAPH_ENABLE_SHARED_MEMORY), 1)
    # Boost.Interprocess is header-only, but its named shared memory objects
    # need librt on Linux (not on macOS, where the equivalent is in libc).
    ifneq ($(shell uname -s), Darwin)
        LIBS += -lrt
    endif
endif

# Apple Clang does not support OpenMP directly, so we need special handling.
ifeq ($(shell uname -s), Darwin)
    ifeq ($(GBWTGRAPH_USE_OPENMP), 1)
    # The compiler complains about -fopenmp instead of missing input.
    ifeq ($(strip $(shell $(MY_CXX) -fopenmp /dev/null -o/dev/null 2>&1 | grep fopenmp | wc -l)), 1)
        $(info The compiler is Apple Clang that needs libomp for OpenMP support.)

        # The compiler only needs to do the preprocessing.
        PARALLEL_FLAGS=-Xpreprocessor -fopenmp -pthread -DGBWTGRAPH_USE_OPENMP -DGBWT_USE_OPENMP

        # Find libomp installed by Homebrew or MacPorts.
        ifeq ($(shell if [ -e $(HOMEBREW_PREFIX)/include/omp.h ]; then echo 1; else echo 0; fi), 1)
            $(info Found libomp installed by Homebrew and linked to $(HOMEBREW_PREFIX).)
            INCLUDES += -I$(HOMEBREW_PREFIX)/include
            LIBS += -L$(HOMEBREW_PREFIX)/lib
        else ifeq ($(shell if [ -d $(HOMEBREW_PREFIX)/opt/libomp/include ]; then echo 1; else echo 0; fi), 1)
            $(info Found a keg-only libomp installed by Homebrew at $(HOMEBREW_PREFIX)/opt/libomp.)
            INCLUDES += -I$(HOMEBREW_PREFIX)/opt/libomp/include
            LIBS += -L$(HOMEBREW_PREFIX)/opt/libomp/lib
        else ifeq ($(shell if [ -d /opt/local/lib/libomp ]; then echo 1; else echo 0; fi), 1)
            $(info Found libomp installed by MacPorts at /opt/local.)
            INCLUDES += -I/opt/local/include/libomp
            LIBS += -L/opt/local/lib/libomp
        else
            $(error Could not find libomp. Please install it using Homebrew or MacPorts)
        endif

        # We also need to link it.
        LIBS += -lomp
    endif
    endif
endif

CXX_FLAGS=$(MY_CXX_FLAGS) $(PARALLEL_FLAGS) $(DEFINES) $(MY_CXX_OPT_FLAGS) $(INCLUDES)

HEADERS=$(wildcard include/gbwtgraph/*.h)
LIBOBJS=$(addprefix $(BUILD_OBJ)/,algorithms.o cached_gbwtgraph.o gbwtgraph.o gbz.o gfa.o index.o internal.o minimizer.o naive_graph.o path_cover.o subgraph.o utils.o)
LIBRARY=$(BUILD_LIB)/libgbwtgraph.a

PROGRAMS=$(addprefix $(BUILD_BIN)/,canonical_gfa chunk_gbz gfa2gbwt gbz_extract gbz_find gbz_stats kmer_freq subgraph_query)
OBSOLETE=gfa2gbwt

.PHONY: all clean directories test
all: directories $(LIBRARY) $(PROGRAMS)

directories: $(BUILD_BIN) $(BUILD_LIB) $(BUILD_OBJ)

$(BUILD_BIN):
	mkdir -p $@

$(BUILD_LIB):
	mkdir -p $@

$(BUILD_OBJ):
	mkdir -p $@

$(BUILD_OBJ)/%.o:$(SOURCE_DIR)/%.cpp $(HEADERS)
	$(MY_CXX) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -c -o $@ $<

$(LIBRARY):$(LIBOBJS)
	ar rcs $@ $(LIBOBJS)

$(BUILD_BIN)/%:$(BUILD_OBJ)/%.o $(LIBRARY)
	$(MY_CXX) $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

test:$(LIBRARY)
	cd tests && $(MAKE) test

clean:
	rm -rf $(BUILD_BIN) $(BUILD_LIB) $(BUILD_OBJ)
	rm -f *.o *.a $(OBSOLETE)
	cd tests && $(MAKE) clean
