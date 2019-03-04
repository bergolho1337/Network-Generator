CXX = g++

all: generate

generate: generate.cpp
	$(CXX) -o $@ $< 

clean:
	rm -f generate
	rm -f files/*

.PHONY: clean
