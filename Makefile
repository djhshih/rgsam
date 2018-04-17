CXX = g++
CPPFLAGS = -O3
DESTDIR ?= /usr/local

all: bin/rgsam
	

bin/rgsam: rgsam.cpp
	mkdir -p bin
	$(CXX) $(CPPFLAGS) $? -o $@

check: rgsam.cpp
	mkdir -p tmp
	$(CXX) -coverage -O0 $? -o tmp/check
	! tmp/check
	! tmp/check fly
	! tmp/check tag -i data/illumina-1.8.sam
	# test qnames
	tmp/check qnames
	# test collect on sam files
	tmp/check collect -q illumina-1.8 -i data/illumina-1.8.sam -s sample1 -l library1 -o tmp/illumina-1.8.sam.rg.txt
	diff data/ans/illumina-1.8.sam.rg.txt tmp/illumina-1.8.sam.rg.txt
	# test collect on fastq files
	tmp/check collect -q illumina-1.0 -i data/illumina-1.0.fq -s sample1 -l library1 -o tmp/illumina-1.0.fq.rg.txt
	diff data/ans/illumina-1.0.fq.rg.txt tmp/illumina-1.0.fq.rg.txt
	tmp/check collect -i data/illumina-1.8.fq -s sample1 -l library1 -o tmp/illumina-1.8.fq.rg.txt
	diff data/ans/illumina-1.8.fq.rg.txt tmp/illumina-1.8.fq.rg.txt
	tmp/check collect -i data/illumina-1.8.fq -s sample1 -l library1 > tmp/illumina-1.8.fq.rg.txt
	diff data/ans/illumina-1.8.fq.rg.txt tmp/illumina-1.8.fq.rg.txt
	cat data/illumina-1.8.fq | tmp/check collect -f fastq -s sample1 -l library1 > tmp/illumina-1.8.fq.rg.txt
	diff data/ans/illumina-1.8.fq.rg.txt tmp/illumina-1.8.fq.rg.txt
	cat data/illumina-1.8.fq | tmp/check collect -i - -o - -f fastq -s sample1 -l library1 > tmp/illumina-1.8.fq.rg.txt
	diff data/ans/illumina-1.8.fq.rg.txt tmp/illumina-1.8.fq.rg.txt
	tmp/check collect -q broad-1.0 -i data/broad-1.0.fq -s sample1 -l library1 -o tmp/broad-1.0.fq.rg.txt
	diff data/ans/broad-1.0.fq.rg.txt tmp/broad-1.0.fq.rg.txt
	# test tag
	tmp/check tag -q illumina-1.8 -i data/illumina-1.8.sam -r data/ans/illumina-1.8.sam.rg.txt -o tmp/illumina-1.8.rg.sam
	diff data/ans/illumina-1.8.rg.sam tmp/illumina-1.8.rg.sam

coverage: check
	gcov rgsam.cpp

test: check
	

install: bin/rgsam
	mkdir -p $(DESTDIR)/bin/
	install bin/rgsam $(DESTDIR)/bin/

clean:
	rm -f bin/rgsam
	rm -f *.exe *.gcov *.gcno *.gcda
	rm -rf tmp

