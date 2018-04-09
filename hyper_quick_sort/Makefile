all: hqs

hqs: hyper_qsort main serial_qsort
	mpicc hyper_qsort.o main.o serial_qsort.o -g -o hqs

hyper_qsort:
	mpicc hyper_qsort.c -c -g -o hyper_qsort.o

main: main.c hyper_qsort
	mpicc -c main.c -g -o main.o

serial_qsort:
	mpicc -c serial_qsort.c -g -o serial_qsort.o
