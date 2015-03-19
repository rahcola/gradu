CFLAGS=-Wall -Werror -pedantic -O3

de-bruijn: src/de-bruijn.cpp
	g++ -std=c++11 $(CFLAGS) -I /usr/local/include -L /usr/local/lib/ $< -o $@ -ldivsufsort -ldivsufsort64 -lsdsl
