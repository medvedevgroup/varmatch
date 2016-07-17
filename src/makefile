CXX=g++
CXXFLAGS=-std=c++11 -pthread -g
CXXFLAGS2=-I ../include
CXXFLAGTBB=-ltbb

all: vm-core

vm-core: vm.cpp decomposematch.cpp util.cpp
	$(CXX) $(CXXFLAGS) $(CXXFLAGS2) -o $@ $^
	cp $@ ../$@

filter_hc: filter_hc.cpp util.cpp
	$(CXX) $(CXXFLAGS) $(CXXFLAGS2) -o $@ $^

filter_cv: filter_cv.cpp util.cpp
	$(CXX) $(CXXFLAGS) $(CXXFLAGS2) -o $@ $^

clean:
	rm -f vm-core
	rm -f *.o
