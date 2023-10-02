SDSL_DIR=../../sdsl-lite
include $(SDSL_DIR)/Make.helper

MAIN_DIR=..
LIBRARY=$(MAIN_DIR)/lib/libgbwtgraph.a

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
            ifeq ($(shell if [ -d $(HOMEBREW_PREFIX)/opt/libomp/include ]; then echo 1; else echo 0; fi), 1)
                # libomp moved to these directories, recently, because it is now keg-only to not fight GCC
                PARALLEL_FLAGS += -I$(HOMEBREW_PREFIX)/opt/libomp/include
                LIBS += -L$(HOMEBREW_PREFIX)/opt/libomp/lib
            else
                PARALLEL_FLAGS += -I$(HOMEBREW_PREFIX)/include
                LIBS += -L$(HOMEBREW_PREFIX)/lib
            endif
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

HEADERS=$(wildcard $(GBWT_DIR)/include/gbwt/*.h) shared.h
PROGRAMS=test_utils test_gbwtgraph test_cached_gbwtgraph test_gfa test_gbz test_minimizer test_index test_algorithms test_path_cover

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
