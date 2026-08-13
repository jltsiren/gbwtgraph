SDSL_DIR=../sdsl-lite
include $(SDSL_DIR)/Make.helper

BUILD_BIN=bin
BUILD_LIB=lib
BUILD_OBJ=obj
SOURCE_DIR=src

# --- Optional build-time behaviors, overridable on the command line (e.g.
# `make GBWTGRAPH_USE_EXCEPTIONS=0`). Defaults match GBWTGraph's traditional
# behavior. See include/gbwtgraph/error_handling.h and
# include/gbwtgraph/gbwtgraph.h / gbz.h (GBZ, GBWTGraph) for what each one
# changes.

# Report fatal errors by throwing C++ exceptions (1) or not (0). Building
# with 0 is needed to embed GBWTGraph in exceptions-free code, such as
# Google's DeepVariant, which builds with Bazel and -fno-exceptions.
GBWTGRAPH_USE_EXCEPTIONS=1
# When GBWTGRAPH_USE_EXCEPTIONS=0, report fatal errors via Abseil's
# ABSL_LOG(FATAL) (1) instead of stderr + std::abort() (0). No effect when
# GBWTGRAPH_USE_EXCEPTIONS=1.
GBWTGRAPH_USE_ABSEIL_LOGGING=0
# Build with OpenMP parallelism (1) or without it (0).
GBWTGRAPH_USE_OPENMP=1
# Enable GBZ/GBWTGraph's shared-memory-backed node sequence storage, which
# needs Boost.Interprocess (see include/gbwtgraph/gbwtgraph.h and gbz.h) and
# a gbwt library that was itself built with GBWT_ENABLE_SHARED_MEMORY=1
# (for gbwt::SharedMemCharAllocatorType). Defaults to on if a Boost
# installation is found, off otherwise; set explicitly to override the
# autodetection either way.
GBWTGRAPH_ENABLE_SHARED_MEMORY=$(shell echo '\#include <boost/interprocess/managed_shared_memory.hpp>' | $(MY_CXX) -std=c++17 -E -x c++ - >/dev/null 2>&1 && echo 1 || echo 0)

# NOTE: gbwt's own headers (transitively included everywhere) make the same
# behavioral choices under their own GBWT_* macros of the same shape, and
# gbwt's declarations must match however lib/libgbwt.a was actually built
# (mismatched macros between a header and the .a it links against produce
# either hard-to-read compile errors, e.g. conflicting `#pragma omp`-guarded
# declarations of omp_get_max_threads() against the real <omp.h>, or linker
# errors). So each GBWTGRAPH_* choice below also sets the matching GBWT_*
# macro, on the assumption -- which the build of this library requires the
# caller to uphold -- that ../gbwt was itself built with matching flags
# (e.g. `make GBWT_USE_EXCEPTIONS=0` alongside `make GBWTGRAPH_USE_EXCEPTIONS=0`).
DEFINES=
ifeq ($(GBWTGRAPH_USE_EXCEPTIONS), 0)
    DEFINES += -DGBWTGRAPH_NO_EXCEPTIONS -DGBWT_NO_EXCEPTIONS
    ifeq ($(GBWTGRAPH_USE_ABSEIL_LOGGING), 1)
        DEFINES += -DGBWTGRAPH_USE_ABSEIL_LOGGING -DGBWT_USE_ABSEIL_LOGGING
    endif
endif

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

ifeq ($(GBWTGRAPH_USE_ABSEIL_LOGGING), 1)
    ifeq ($(shell pkg-config --exists absl_log absl_log_internal_message && echo 1), 1)
        $(info Found Abseil logging.)
        INCLUDES += $(shell pkg-config --cflags absl_log absl_log_internal_message)
        LIBS += $(shell pkg-config --libs absl_log absl_log_internal_message)
    else
        $(error GBWTGRAPH_USE_ABSEIL_LOGGING=1 but could not find Abseil logging via pkg-config. Please update PKG_CONFIG_PATH)
    endif
endif

ifeq ($(GBWTGRAPH_ENABLE_SHARED_MEMORY), 1)
    DEFINES += -DGBWTGRAPH_ENABLE_SHARED_MEMORY
    # gbwt's own header declarations (e.g. StringArray's shared-memory
    # constructors, and the SharedMemCharAllocatorType typedef we reuse for
    # GBZ/GBWTGraph) are gated by GBWT_ENABLE_SHARED_MEMORY, which must be
    # defined here to match however lib/../gbwt/lib/libgbwt.a itself was
    # actually built -- otherwise the two sides disagree about StringArray's
    # constructor signatures and gbwtgraph fails to compile or link.
    DEFINES += -DGBWT_ENABLE_SHARED_MEMORY
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
