
# Makefile for Branching Process Simulation Project
# homebrew c++ compiler (in /usr/local/Cellar/gcc/4.9.1)
CC=g++

# compiler flags
CFLAGS=-W -Wno-unused-parameter -Wno-long-long -pedantic -Wno-variadic-macros -std=c++11

INCLUDE=

all: basic_gillespie.o target

test: basic_gillespie.o testing

basic_gillespie.o: basic_gillespie.cpp basic_gillespie.h timekeeper.h
	$(CC) $(CFLAGS) $(INCLUDE) -c basic_gillespie.cpp

target: yan_network_cat.cpp basic_gillespie.o
	$(CC) $(CFLAGS) $(INCLUDE) yan_network_cat.cpp basic_gillespie.o -o yanc

testing: testing.cpp basic_gillespie.o
	$(CC) $(CFLAGS) $(INCLUDE) testing.cpp basic_gillespie.o -o testing

sub: yan_subnetwork.cpp basic_gillespie.o
	$(CC) $(CFLAGS) $(INCLUDE) yan_subnetwork.cpp basic_gillespie.o -o ysub

clean: 
	rm -rf *.o testing yan *~
