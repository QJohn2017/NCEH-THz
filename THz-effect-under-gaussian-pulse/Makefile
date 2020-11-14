CPP_FLAGS = -std=c++17 -O3
all:
	g++ main_o_THz.cc -O3 -o oTHz -std=c++17
	g++ main_w_THz.cc -O3 -o wTHz -std=c++17
	g++ main_pulse_THz.cc -O3 -o pulse_THz -std=c++17
	g++ main_pulse_gaussian.cc -O3 -o pulse_gaussian -std=c++17
	g++ ../erf.cc -c -o cerf.o $(CPP_FLAGS)
	g++ main_w_THz_analytic.cc -c -o wTHz_ana.o $(CPP_FLAGS)
	g++ cerf.o wTHz_ana.o -O3 -o wTHz_ana
clean:
	rm -f *~
