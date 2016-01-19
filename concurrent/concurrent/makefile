CXX=g++
CXXFLAGS=-std=c++11 -pthread

all: concurrent

concurrent: concurrent.cpp vcf.cpp util.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -f concurrent
