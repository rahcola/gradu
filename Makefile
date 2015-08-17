SHELL=/bin/sh
CXX=g++
CXXFLAGS=-c -std=c++11 -Wall -Werror -pedantic -O3
LDFLAGS=-ldivsufsort -ldivsufsort64 -lsdsl

all: target/de-bruijn

target/%.o: src/%.cpp | target
	$(CXX) -c $(CXXFLAGS) $< -o $@

target/de-bruijn: target/de-bruijn.o target/BidirectionalBWTIndex.o target/DeBruijn.o | target
	$(CXX) -o $@ $^ $(LDFLAGS)

target:
	mkdir target

clean:
	rm -f ./target/*
