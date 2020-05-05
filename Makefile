# Use C++11, dont warn on long-to-float conversion
CXXFLAGS += -std=c++11 -Wno-conversion

# Default to using system's default version of python
PYTHON_BIN     ?= python
PYTHON_CONFIG  := $(PYTHON_BIN)-config
PYTHON_INCLUDE ?= $(shell $(PYTHON_CONFIG) --includes)
EXTRA_FLAGS    := $(PYTHON_INCLUDE)
# NOTE: Since python3.8, the correct invocation is `python3-config --libs --embed`. 
# So of course the proper way to get python libs for embedding now is to
# invoke that, check if it crashes, and fall back to just `--libs` if it does.
LDFLAGS        += $(shell if $(PYTHON_CONFIG) --libs --embed >/dev/null; then $(PYTHON_CONFIG) --libs --embed; else $(PYTHON_CONFIG) --libs; fi)

# Either finds numpy or set -DWITHOUT_NUMPY
#EXTRA_FLAGS     += $(shell $(PYTHON_BIN) $(CURDIR)/numpy_flags.py)
WITHOUT_NUMPY   := $(findstring $(EXTRA_FLAGS), WITHOUT_NUMPY)

# Examples requiring numpy support to compile
EXAMPLES        := basicDM basicDMverbose basicNRMverbose basicNRM predatorpreyDM  predatorpreyDMverbose predatorpreyNRM predatorpreyNRMverbose

# Prefix every example with 'examples/build/'
EXAMPLE_TARGETS := $(patsubst %,./build/%,$(EXAMPLES))

.PHONY: examples

examples: $(EXAMPLE_TARGETS)

docs:
	doxygen
	moxygen doc/xml --noindex -o doc/api.md

# Assume every *.cpp file is a separate example
$(EXAMPLE_TARGETS): ./build/%: ./examples/%.cpp matplotlibcpp.h crn.h crn.cpp
	mkdir -p ./build
	$(CXX) -o $@ $<  crn.cpp $(EXTRA_FLAGS) $(CXXFLAGS) $(LDFLAGS)

clean:
	rm -f ${EXAMPLE_TARGETS}
