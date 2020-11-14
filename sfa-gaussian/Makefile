CPP_FLAGS = -std=c++17 -O3

all:
	g++ ../erf.cc -c -o cerf.o $(CPP_FLAGS)
	g++ ./sfa.cc -c -o sfa.o $(CPP_FLAGS)
	g++ cerf.o sfa.o -O3 -o app
clean:
	rm -f *~
