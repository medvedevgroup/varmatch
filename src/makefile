CXX=g++
CXXFLAGS=-std=c++11 -pthread -g
CXXFLAGS2=-I ../include

all: vm-core vm-convert

vm-convert: formatconverter.cpp util.cpp
	$(CXX) $(CXXFLAGS) $(CXXFLAGS2) -o $@ $^

vm-core: vm_giab.cpp mergevariants.cpp util.cpp
	$(CXX) $(CXXFLAGS) $(CXXFLAGS2) -o $@ $^
#	cp $@ ../$@

filter_hc: filter_hc.cpp util.cpp
	$(CXX) $(CXXFLAGS) $(CXXFLAGS2) -o $@ $^

filter_cv: filter_cv.cpp util.cpp
	$(CXX) $(CXXFLAGS) $(CXXFLAGS2) -o $@ $^

clean:
	rm -f vm-core
	rm -f *.o
	rm -f vm-convert
