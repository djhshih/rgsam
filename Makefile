CXX=g++
CPPFLAGS=-O3

all: bin/rgsam
	

bin/rgsam: rgsam.cpp
	mkdir -p bin
	$(CXX) $(CPPFLAGS) $? -o $@

check: rgsam.cpp
	mkdir -p tmp
	$(CXX) -coverage -O0 $? -o tmp/check
	! tmp/check
	! tmp/check fly
	# test collect
	tmp/check collect data/illumina-1.8.sam sample1 library1 tmp/illumina-1.8.sam.rg.txt
	diff data/ans/illumina-1.8.sam.rg.txt tmp/illumina-1.8.sam.rg.txt
	# test collectfq
	tmp/check collectfq data/illumina-1.8.fq sample1 library1 tmp/illumina-1.8.fq.rg.txt
	diff data/ans/illumina-1.8.fq.rg.txt tmp/illumina-1.8.fq.rg.txt
	# test tag
	tmp/check tag data/illumina-1.8.sam data/ans/illumina-1.8.sam.rg.txt tmp/illumina-1.8.rg.sam
	diff data/ans/illumina-1.8.rg.sam tmp/illumina-1.8.rg.sam

coverage: check
	gcov rgsam.cpp

test: check
	

clean:
	rm -f bin/rgsam
	rm -f *.exe *.gcov *.gcno *.gcda
	rm -rf tmp

