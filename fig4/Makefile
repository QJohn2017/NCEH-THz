CPP_FLAGS = -std=c++17 -O3

all:
	g++ ../erf.cc -c -o cerf.o $(CPP_FLAGS)
	g++ ./main.cc -c -o main.o $(CPP_FLAGS)
	g++ cerf.o main.o -O3 -o app
clean:
	rm -f *~
