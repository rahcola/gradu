CC=g++
CFLAGS=-c -std=c++11 -Wall -Werror -pedantic -O3
LFLAGS=-ldivsufsort -ldivsufsort64 -lsdsl

de-bruijn: de-bruijn.o BidirectionalBWTIndex.o DeBruijn.o
	$(CC) $^ -o $@ $(LFLAGS)

de-bruijn.o: src/de-bruijn.cpp src/BidirectionalBWTIndex.hpp src/DeBruijn.hpp
	$(CC) $(CFLAGS) $< -o $@

BidirectionalBWTIndex.o: src/BidirectionalBWTIndex.cpp src/BidirectionalBWTIndex.hpp
	$(CC) $(CFLAGS) $< -o $@

DeBruijn.o: src/DeBruijn.cpp src/DeBruijn.hpp src/BidirectionalBWTIndex.hpp
	$(CC) $(CFLAGS) $< -o $@
