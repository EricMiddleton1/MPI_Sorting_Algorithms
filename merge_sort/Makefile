all: sort

sort: merge_sort main serial_qsort
	mpicc merge_sort.o main.o serial_qsort.o -g -o sort

merge_sort:
	mpicc merge_sort.c -c -g -o merge_sort.o

main: main.c merge_sort
	mpicc -c main.c -g -o main.o

serial_qsort:
	mpicc -c serial_qsort.c -g -o serial_qsort.o
