SDSL_DIR=../../sdsl-lite
include $(SDSL_DIR)/Make.helper

MAIN_DIR=..
LIBRARY=$(MAIN_DIR)/libgbwtgraph.a

# Multithreading with OpenMP.
PARALLEL_FLAGS=-fopenmp -pthread

LIBS=-L$(LIB_DIR) -lgbwt -lhandlegraph -lsdsl -ldivsufsort -ldivsufsort64 -lgtest -lgtest_main

# Apple Clang does not support OpenMP directly, so we need special handling.
ifeq ($(shell uname -s), Darwin)
    # The compiler complains about -fopenmp instead of missing input.
    ifeq ($(strip $(shell $(MY_CXX) -fopenmp /dev/null -o/dev/null 2>&1 | grep fopenmp | wc -l)), 1)
        # The compiler only needs to do the preprocessing.
        PARALLEL_FLAGS = -Xpreprocessor -fopenmp -pthread

        # If HOMEBREW_PREFIX is specified, libomp probably cannot be found automatically.
        ifdef HOMEBREW_PREFIX
            PARALLEL_FLAGS += -I$(HOMEBREW_PREFIX)/include
            LIBS += -L$(HOMEBREW_PREFIX)/lib
        # Macports installs libomp to /opt/local/lib/libomp
        else ifeq ($(shell if [ -d /opt/local/lib/libomp ]; then echo 1; else echo 0; fi), 1)
            PARALLEL_FLAGS += -I/opt/local/include/libomp
            LIBS += -L/opt/local/lib/libomp
        endif

        # We also need to link it.
        LIBS += -lomp
    endif
endif

CXX_FLAGS=$(MY_CXX_FLAGS) $(PARALLEL_FLAGS) $(MY_CXX_OPT_FLAGS) -I$(MAIN_DIR)/include -I$(INC_DIR)
SOURCES=$(wildcard *.cpp)
HEADERS=$(wildcard ../include/gbwtgraph/*.h)
HEADERS+=shared.h
OBJS=$(SOURCES:.cpp=.o)
PROGRAMS=test_algorithms test_cached_gbwtgraph test_gbwtgraph test_gfa test_index test_minimizer test_path_cover test_utils

all:$(PROGRAMS)

%.o:%.cpp $(HEADERS)
	$(MY_CXX) $(CXX_FLAGS) -c $<

test_algorithms:test_algorithms.o $(LIBRARY) $(HEADERS)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

test_cached_gbwtgraph:test_cached_gbwtgraph.o $(LIBRARY) $(HEADERS)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

test_gbwtgraph:test_gbwtgraph.o $(LIBRARY) $(HEADERS)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

test_gfa:test_gfa.o $(LIBRARY) $(HEADERS)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

test_index:test_index.o $(LIBRARY) $(HEADERS)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

test_minimizer:test_minimizer.o $(LIBRARY) $(HEADERS)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

test_path_cover:test_path_cover.o $(LIBRARY) $(HEADERS)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

test_utils:test_utils.o $(LIBRARY) $(HEADERS)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

test:$(PROGRAMS)
	$(foreach PROG,$(PROGRAMS),./$(PROG) || exit 1;)

clean:
	rm -f $(PROGRAMS) $(OBJS)
