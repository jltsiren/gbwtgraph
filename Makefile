SDSL_DIR=../sdsl-lite
include $(SDSL_DIR)/Make.helper

# Multithreading with OpenMP.
PARALLEL_FLAGS=-fopenmp -pthread

LIBS=-L$(LIB_DIR) -lgbwt -lhandlegraph -lsdsl -ldivsufsort -ldivsufsort64

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

CXX_FLAGS=$(MY_CXX_FLAGS) $(PARALLEL_FLAGS) $(MY_CXX_OPT_FLAGS) -Iinclude -I$(INC_DIR)
LIBOBJS=algorithms.o cached_gbwtgraph.o gbwtgraph.o gfa.o internal.o minimizer.o path_cover.o utils.o
SOURCES=$(wildcard *.cpp)
HEADERS=$(wildcard include/gbwtgraph/*.h)
OBJS=$(SOURCES:.cpp=.o)

LIBRARY=libgbwtgraph.a
PROGRAMS=gfa2gbwt

all:$(LIBRARY) $(PROGRAMS)

%.o:%.cpp $(HEADERS)
	$(MY_CXX) $(CXX_FLAGS) -c $<

gfa2gbwt:gfa2gbwt.o $(LIBRARY)
	$(MY_CXX) $(LDFLAGS) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

$(LIBRARY):$(LIBOBJS)
	ar rcs $@ $(LIBOBJS)

test:$(LIBRARY)
	cd tests && $(MAKE) test

clean:
	rm -f $(OBJS) $(LIBRARY) $(PROGRAMS)
