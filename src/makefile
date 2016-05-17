CXX=g++
CXXFLAGS=-std=c++11 -pthread -g
CXXFLAGS2=-I ./include
CXXFLAGTBB=-ltbb

all: vm

vm: vm.cpp vcf.cpp removeduplicate.cpp diploid.cpp wholegenome.cpp util.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^
	cp $@ ../$@

split: splitvcf.cpp util.cpp
	$(CXX) $(CXXFLAGS) $(CXXFLAGS2) -o $@ $^

clean:
	rm -f vm
	rm -f *.o
	rm -f split
