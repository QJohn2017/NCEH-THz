all:
	cd ../sfa-mono-THz/ && make && ./wTHz # generate para_1.dat data_1.dat field_1.dat
	cd ../ana-mono-THz/ && make && ./app  # generate res.dat
	cd ../ana-mono && make && ./app       # generate res.dat
	python plot.py