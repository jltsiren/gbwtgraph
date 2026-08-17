SDSL_DIR=../../sdsl-lite
include $(SDSL_DIR)/Make.helper

MAIN_DIR=..
LIBRARY=$(MAIN_DIR)/lib/libgbwtgraph.a

# This must be set the same way as when ../lib/libgbwtgraph.a itself was
# built (e.g. `make GBWTGRAPH_USE_OPENMP=0` in both places), since it
# selects which declarations and macro expansions the tests see in
# gbwtgraph's headers. See ../Makefile for what it does.
GBWTGRAPH_USE_OPENMP=1
# Not user-facing like the flag above: see GBWTGRAPH_ENABLE_SHARED_MEMORY
# in ../Makefile for why this is derived from gbwt/config.h instead.
GBWTGRAPH_ENABLE_SHARED_MEMORY=$(shell echo '\#include <gbwt/config.h>' | $(MY_CXX) -I$(INC_DIR) -E -x c++ - 2>/dev/null | grep -q GBWT_ENABLE_SHARED_MEMORY && echo 1 || echo 0)

DEFINES=

# Multithreading with OpenMP.
ifeq ($(GBWTGRAPH_USE_OPENMP), 1)
    PARALLEL_FLAGS=-fopenmp -pthread -DGBWTGRAPH_USE_OPENMP -DGBWT_USE_OPENMP
else
    PARALLEL_FLAGS=-pthread -Wno-unknown-pragmas
endif

# Directories for dependencies.
INCLUDES=-I$(MAIN_DIR)/include -I$(INC_DIR)
LIBS=-L$(LIB_DIR) -lgbwt -lhandlegraph -lsdsl -lgtest -lgtest_main

# Use pkg-config to find system dependencies.
ifeq ($(shell pkg-config --exists libcrypto --exists libzstd && echo 1), 1)
    $(info Found libcrypto and libzstd.)
    INCLUDES += $(shell pkg-config --cflags libcrypto libzstd)
    LIBS += $(shell pkg-config --libs libcrypto libzstd)
else
    $(error Could not find libcrypto or libzstd. Please update PKG_CONFIG_PATH)
endif

ifeq ($(GBWTGRAPH_ENABLE_SHARED_MEMORY), 1)
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

    # gtest isn't in pkg-config but also we can't find its libs by default if
    # they are installed via Homebrew.
    ifeq ($(shell if [ -e $(HOMEBREW_PREFIX)/lib/libgtest.a ]; then echo 1; else echo 0; fi), 1)
        $(info Found gtest installed by Homebrew and linked to $(HOMEBREW_PREFIX).)
        INCLUDES += -I$(HOMEBREW_PREFIX)/include
        LIBS += -L$(HOMEBREW_PREFIX)/lib
    endif
endif

CXX_FLAGS=$(MY_CXX_FLAGS) $(PARALLEL_FLAGS) $(DEFINES) $(MY_CXX_OPT_FLAGS) $(INCLUDES)

HEADERS=$(wildcard $(MAIN_DIR)/include/gbwtgraph/*.h) shared.h
PROGRAMS=test_utils test_gbwtgraph test_cached_gbwtgraph test_gfa test_gbz test_minimizer test_index test_algorithms test_path_cover test_subgraph test_naive_graph

.PHONY: all clean test
all:$(PROGRAMS)

%.o:%.cpp $(HEADERS)
	$(MY_CXX) $(CXX_FLAGS) -c $<

$(PROGRAMS):%:%.o $(LIBRARY)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

test:$(PROGRAMS)
	$(foreach PROG,$(PROGRAMS),./$(PROG) || exit 1;)

clean:
	rm -f $(PROGRAMS) *.o
