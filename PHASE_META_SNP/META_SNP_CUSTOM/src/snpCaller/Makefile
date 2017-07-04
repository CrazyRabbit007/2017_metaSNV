include ../../SETUPFILE
# Using gcc for now
CC=g++
# Flags for the compiles
CFLAGS=-c -O2
# Path to BOOST. This is important.
#BOOST=/g/bork8/costea/boost_1_53_0
#BOOST=/g/bork5/rmuench/gitlab/boost_1_55_0
INCLUDE=.

all: snpCall

snpCall: snpCall.o
	$(CC) snpCall.o -o snpCall

snpCall.o: call_vC.cpp
	$(CC) -I$(BOOST) -I$(INCLUDE) $(CFLAGS) $< -o $@

clean:
	rm -rf *o snpCall
