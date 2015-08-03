CC=g++
CFLAGS=-c -std=c++11 -Wall -Werror -pedantic -O3
LFLAGS=-ldivsufsort -ldivsufsort64 -lsdsl

de-bruijn: de-bruijn.o BidirectionalBWTIndex.o
	$(CC) $^ -o $@ $(LFLAGS)

de-bruijn.o: src/de-bruijn.cpp src/BidirectionalBWTIndex.hpp
	$(CC) $(CFLAGS) $< -o $@

BidirectionalBWTIndex.o: src/BidirectionalBWTIndex.cpp src/BidirectionalBWTIndex.hpp
	$(CC) $(CFLAGS) $< -o $@
