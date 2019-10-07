SDSL_DIR=../../sdsl-lite
include $(SDSL_DIR)/Make.helper

MAIN_DIR=..
LIBRARY=$(MAIN_DIR)/libgbwtgraph.a

# Multithreading with OpenMP.
PARALLEL_FLAGS=-fopenmp -pthread

LIBS=-L$(LIB_DIR) -lgbwt -lhandlegraph -lsdsl -ldivsufsort -ldivsufsort64 -lgtest -lgtest_main

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

CXX_FLAGS=$(MY_CXX_FLAGS) $(PARALLEL_FLAGS) $(MY_CXX_OPT_FLAGS) -I$(MAIN_DIR)/include -I$(INC_DIR)
SOURCES=$(wildcard *.cpp)
HEADERS=$(wildcard ../include/gbwtgraph/*.h)
HEADERS+=shared.h
OBJS=$(SOURCES:.cpp=.o)
PROGRAMS=test_gbwtgraph test_gfa test_index test_minimizer

all:$(PROGRAMS)

%.o:%.cpp $(HEADERS)
	$(MY_CXX) $(CXX_FLAGS) -c $<

test_gbwtgraph:test_gbwtgraph.o $(LIBRARY) $(HEADERS)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

test_gfa:test_gfa.o $(LIBRARY) $(HEADERS)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

test_index:test_index.o $(LIBRARY) $(HEADERS)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

test_minimizer:test_minimizer.o $(LIBRARY) $(HEADERS)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

test:$(PROGRAMS)
	$(foreach PROG,$(PROGRAMS),./$(PROG) || exit 1;)

clean:
	rm -f $(PROGRAMS) $(OBJS)
