SDSL_DIR=../sdsl-lite
include $(SDSL_DIR)/Make.helper

# Multithreading with OpenMP.
PARALLEL_FLAGS=-fopenmp -pthread

LIBS=-L$(LIB_DIR) -lgbwt -lhandlegraph -lsdsl -ldivsufsort -ldivsufsort64

ifeq ($(shell uname -s),Darwin)
    # Our compiler might be clang that lacks -fopenmp support.
    # Sniff that
    ifeq ($(strip $(shell $(MY_CXX) -fopenmp /dev/null -o/dev/null 2>&1 | grep fopenmp | wc -l)), 1)
		# The compiler complained about fopenmp instead of its nonsense input file.
        # We need to use the hard way of getting OpenMP not bundled with the compiler.
        
        # The compiler only needs to do the preprocessing
        PARALLEL_FLAGS=-Xpreprocessor -fopenmp -pthread

        ifeq ($(shell if [ -d /opt/local/lib/libomp ];then echo 1;else echo 0;fi), 1)
            # Use /opt/local/lib/libomp if present, because Macports installs libomp there.
            # Brew is supposed to put it somewhere the compiler can find it by default.
            LIBS += -L/opt/local/lib/libomp
            # And we need to find the includes. Homebrew puts them in the normal place
            # but Macports hides them in "libomp"
            PARALLEL_FLAGS += -I/opt/local/include/libomp
        endif

        # We also need to link it
        LIBS += -lomp
    endif
endif

CXX_FLAGS=$(MY_CXX_FLAGS) $(PARALLEL_FLAGS) $(MY_CXX_OPT_FLAGS) -Iinclude -I$(INC_DIR)
LIBOBJS=algorithms.o cached_gbwtgraph.o gbwtgraph.o gfa.o minimizer.o path_cover.o utils.o
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
