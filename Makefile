CXX = g++
CXXFLAGS = -g -O2 -std=c++11
LDFLAGS =

.PHONY: clean

all: bin/bin_coverage

bin/bin_coverage: src/bin_coverage.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f bin/bin_coverage