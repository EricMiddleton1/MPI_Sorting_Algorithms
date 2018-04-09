all: sort

sort: psrs main serial_qsort
	mpicc psrs.o main.o serial_qsort.o -g -o sort

psrs:
	mpicc psrs.c -c -g -o psrs.o

main: main.c psrs
	mpicc -c main.c -g -o main.o

serial_qsort:
	mpicc -c serial_qsort.c -g -o serial_qsort.o
