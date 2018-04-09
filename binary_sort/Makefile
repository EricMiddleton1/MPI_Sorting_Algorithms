all: sort

sort: binary_sort main serial_binary_sort sort_util
	mpicc binary_sort.o main.o serial_binary_sort.o sort_util.o -g -o sort

binary_sort:
	mpicc binary_sort.c -c -g -o binary_sort.o

main: main.c
	mpicc -c main.c -g -o main.o

serial_binary_sort:
	mpicc -c serial_binary_sort.c -g -o serial_binary_sort.o

sort_util: sort_util.c
	mpicc -c sort_util.c -g -o sort_util.o

clean:
	rm *.o
