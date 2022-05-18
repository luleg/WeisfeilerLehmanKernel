all: main rmo

main: wlstkernel md5
	g++ -O3 -fopenmp -lrt -o wlkernel wlkernel.cpp WeisfeilerLehmanSubtreeKernel.o md5.o -I../utils -I../hash-library-master

wlstkernel:
	g++ -O3 -fopenmp -c ../utils/WeisfeilerLehmanSubtreeKernel.cpp -I../hash-library-master

md5:
	g++ -O3 -I../hash-library-master -c ../hash-library-master/md5.cpp

rmo:
	rm *.o
