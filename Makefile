CPPF=-g -O2 -W -Wall -Wno-unused-parameter -Wpointer-arith -Wshadow -Wundef -std=c++11
FSD=src/fstats/

all: bin bin/fasta_stats bin/assem_dir_to_csv

bin:
	mkdir bin

bin/fasta_stats: ${FSD}fasta_stats.cc ${FSD}itoa.cc ${FSD}open_compressed.h
	g++ ${CPPF} -o bin/fasta_stats ${FSD}*.cc

bin/assem_dir_to_csv: scripts/pull_assembly_info.py
	chmod ugo+x scripts/pull_assembly_info.py; \
	cd bin; \
	ln -s ../scripts/pull_assembly_info.py assem_dir_to_csv;\
	cd ..
