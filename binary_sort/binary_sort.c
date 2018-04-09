#include "binary_sort.h"

#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <mpi.h>

static void binary_sort_rec(int arr[], int arr_start, int arr_end, int my_rank, int p_start,
	int p_end);

static void merge(int arr[], int arr_size, int arr_2[], int arr_2_size);

void binary_sort(int arr[], size_t size, int my_rank, int comm_sz) {
	binary_sort_rec(arr, 0, size, my_rank, 0, comm_sz);
}

void binary_sort_rec(int arr[], int arr_start, int arr_end, int my_rank, int p_start, int p_end) {
	MPI_Status status;

	if((p_end - p_start) <= 1) {
		//Base case
		//Perform binary sort on local list
		serial_binary_sort(arr + arr_start, arr_end - arr_start);
	}
	else  {
		int split = p_start + (p_end - p_start)/2 + ((p_end - p_start) % 2);
		int lower_p_start = p_start, lower_p_end = split,
			upper_p_start = split, upper_p_end = p_end;
		int arr_split = arr_start + (arr_start + arr_end)/2;
		int *scratch = NULL;
		int split_size = arr_end - arr_split;

		//Split array in half and send to other process
		if(my_rank == p_start) {
			//Send upper half of array
			MPI_Send(arr + arr_start + arr_split, arr_end - arr_split, MPI_INT, split,
				0, MPI_COMM_WORLD);
		}
		else if(my_rank == split) {
			scratch = (int*)malloc(split_size * sizeof(int));
			
			//Receive upper half of array
			MPI_Recv(scratch, split_size, MPI_INT, p_start, 0, MPI_COMM_WORLD, &status);
		}

		//Recurse
		if(my_rank < split) {
			//Lower-half processes
			binary_sort_rec(arr, arr_start, arr_split, my_rank, p_start, split);
		}
		else {
			//Upper-half processes
			binary_sort_rec(scratch, 0, split_size, my_rank, split, p_end);
		}

		//Merge both sorted halves
		if(my_rank == p_start) {
			scratch = (int*)malloc(split_size * sizeof(int));

			//Receive sorted half
			MPI_Recv(scratch, split_size, MPI_INT, split, 0, MPI_COMM_WORLD, &status);

			//Merge both sorted arrays
			merge(arr + arr_start, arr_split - arr_start, scratch, split_size);

			free(scratch);
		}
		else if(my_rank == split) {
			//Send sorted half
			MPI_Send(scratch, split_size, MPI_INT, p_start, 0, MPI_COMM_WORLD);

			free(scratch);
		}
	}
}


void merge(int arr[], int arr_size, int arr_2[], int arr_2_size) {
	int *scratch = (int*)malloc((arr_size + arr_2_size) * sizeof(int));

	int i_scratch = 0, i_arr = 0, i_arr_2 = 0;

	for(; (i_arr < arr_size) && (i_arr_2 < arr_2_size); ++i_scratch) {
		if(arr[i_arr] < arr_2[i_arr_2]) {
			scratch[i_scratch] = arr[i_arr++];
		}
		else {
			scratch[i_scratch] = arr_2[i_arr_2++];
		}
	}
	if(i_arr < arr_size) {
		memcpy(scratch + i_scratch, arr + i_arr, (arr_size - i_arr) * sizeof(int));
	}
	else if(i_arr_2 < arr_2_size) {
		memcpy(scratch + i_scratch, arr_2 + i_arr_2, (arr_2_size - i_arr_2) * sizeof(int));
	}

	memcpy(arr, scratch, (arr_size + arr_2_size) * sizeof(int));

	free(scratch);
}
