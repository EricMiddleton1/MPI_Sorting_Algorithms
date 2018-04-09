#include <stdio.h>
#include <stdlib.h>
#include <mpi.h> 

#include "binary_sort.h"
#include "sort_util.h"

#define ARRAY_SIZE		1024

void print_array(int* arr, size_t size);

int main(void) {
   int my_rank, comm_sz;

   MPI_Init(NULL, NULL); 
   MPI_Comm_size(MPI_COMM_WORLD, &comm_sz); 
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 

	int* arr;
	if(my_rank == 0) {
		arr = (int*)malloc(ARRAY_SIZE * sizeof(int));
		int i;
		for(i = 0; i < ARRAY_SIZE; ++i) {
			arr[i] = rand() % 100;
		}
	}

	binary_sort(arr, ARRAY_SIZE, my_rank, comm_sz);
	
	if(my_rank == 0) {
		if(validate(arr, ARRAY_SIZE)) {
			printf("[Info] Validation successful!\n");
		}
		else {
			printf("[Error] Validation not successful :(\n");
			print_array(arr, ARRAY_SIZE);
		}

		free(arr);
	}

   MPI_Finalize();
   return 0;
}  /* main */

void print_array(int* arr, size_t size) {
	printf("\t");
	int i;
	for(i = 0; i < size; ++i) {
		printf("%d ", arr[i]);
	}
	printf("\n");
}
