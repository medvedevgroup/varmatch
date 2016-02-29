CXX=g++
CXXFLAGS=-g -std=c++11 -pthread

all: vm

vm: vm.cpp vcf.cpp util.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -f vm
