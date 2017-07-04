include ../../SETUPFILE
# Using gcc for now
CC=g++
# Flags for the compiles
CFLAGS=-c -Wno-sign-compare
INC= -I$(SAMTOOLS) -I$(HTSLIB)
LIBS= -L$(SAMTOOLS) -L$(HTSLIB) -lbam -lhts -pthread -lz

all: qaCompute

qaCompute: qaCompute.o
	$(CC) qaCompute.o -o qaCompute $(LIBS)

qaCompute.o: qaCompute.c
	$(CC) $(INC) $(CFLAGS) -o qaCompute.o qaCompute.c

clean:
	rm -rf *o qaCompute

