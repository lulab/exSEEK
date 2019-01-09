CXX = g++
CXXFLAGS = -g -O2
LDFLAGS =

.PHONY: clean

all: bin/tbed2gbed bin/bin_coverage

bin/tbed2gbed: src/tbed2gbed.o src/formats.o
	$(CXX) -o $@ $^ $(LDFLAGS)

src/tbed2gbed.o: src/tbed2gbed.cpp src/utils.h src/formats.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

src/formats.o: src/formats.cpp src/formats.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

src/bin_coverage.o: src/bin_coverage.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

bin/bin_coverage: src/bin_coverage.o
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	rm -f bin/tbed2gbed bin/bin_coverge
	rm -f src/*.o