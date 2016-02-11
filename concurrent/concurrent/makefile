CXX=g++
CXXFLAGS=-g -std=c++11 -pthread

all: concurrent

concurrent: concurrent.cpp vcf.cpp vcf.h util.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -f concurrent
