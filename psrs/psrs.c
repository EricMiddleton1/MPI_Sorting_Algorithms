#include "psrs.h"
#include "serial_qsort.h"

#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <mpi.h>

static int partition(int arr[], int start, int end, int pivot);
static int merge(int arr[], int *sublists[], int list_counts[], int n_lists);
static int min_index(int *values, int *mask, int n);

void psrs(int arr[], size_t size, int my_rank, int comm_sz) {
	int *my_arr = (int*)malloc(size * sizeof(int)),
		*recv_arr = (int*)malloc(size/comm_sz * sizeof(int)),
		**sublists = (int**)malloc(comm_sz * sizeof(int*)),
		*sub_counts = (int*)malloc(comm_sz * sizeof(int)),
		*pivots = (int*)malloc(comm_sz * sizeof(int));
	int *recv_counts, *displacements;
	int count, i;

	MPI_Request *requests = (MPI_Request*)malloc((comm_sz - 1) * sizeof(MPI_Request));
	MPI_Status status;

	for(i = 0; i < comm_sz; ++i) {
		sublists[i] = (int*)malloc(size/comm_sz * sizeof(int));
	}

	//Distribute partial lists to all processes
	if(my_rank == 0) {
		recv_counts = (int*)malloc(comm_sz * sizeof(int));
		displacements = (int*)malloc(comm_sz * sizeof(int));

		//Send array chunks to other processes
		int i;
		for(i = 1; i < comm_sz; ++i) {
			size_t start = i*size/comm_sz,
				end = (i+1)*size/comm_sz;

			MPI_Send(arr + start, end - start, MPI_INT, i, 0, MPI_COMM_WORLD);
		}

		count = size/comm_sz;
		memcpy(my_arr, arr, count * sizeof(int));
	}
	else {
		//Receive array chunk from master
		my_arr = (int*)malloc(size * sizeof(int));
		MPI_Status status;
		MPI_Recv(my_arr, size, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_INT, &count);
	}

	//Each process sorts partial list
	serial_qsort(my_arr, count);

	//Generate local regular samples
	int* samples = (int*)malloc(comm_sz * sizeof(int));
	int p_sqr = comm_sz * comm_sz;
	for(i = 0; i < comm_sz; ++i) {
		int sample = i*size/p_sqr;
		samples[i] = my_arr[sample];
	}

	//Gather all samples onto root
	int* all_samples;
	if(my_rank == 0) {
		all_samples = (int*)malloc(p_sqr * sizeof(int));
	}
	MPI_Gather(samples, comm_sz, MPI_INT, all_samples, comm_sz, MPI_INT, 0, MPI_COMM_WORLD);

	//Select pivot values
	if(my_rank == 0) {
		for(i = 1; i < comm_sz; ++i) {
			pivots[i-1] = all_samples[i*comm_sz];
		}
	}

	//Broadcast pivot values
	MPI_Bcast(pivots, comm_sz - 1, MPI_INT, 0, MPI_COMM_WORLD);

	//Send sublists
	int list_start = 0;
	for(i = 0; i < comm_sz; ++i) {
		int pivot = (i == (comm_sz - 1)) ? INT_MAX : pivots[i];
		int list_end = partition(my_arr, list_start, count, pivot);

		if(i != my_rank) {
			int i_request = i - (i > my_rank);
			//Non-blocking send
			MPI_Isend(my_arr + list_start, list_end - list_start, MPI_INT, i, 0, MPI_COMM_WORLD,
				&requests[i_request]);
		}
		else {
			//Copy own sublist into array
			sub_counts[i] = list_end - list_start;
			memcpy(sublists[i], my_arr + list_start, sub_counts[i] * sizeof(int));
		}
		
		list_start = list_end;
	}

	//Receive sublists
	for(i = 0; i < comm_sz; ++i) {
		if(i != my_rank) {
			//Blocking receive
			MPI_Recv(sublists[i], size/comm_sz*sizeof(int), MPI_INT, i, 0, MPI_COMM_WORLD, &status);
			MPI_Get_count(&status, MPI_INT, &sub_counts[i]);
		}
	}

	//Wait for non-blocking sends to complete
	MPI_Waitall(comm_sz - 1, requests, MPI_STATUSES_IGNORE);

	//Merge all sublists into sorted list
	count = merge(my_arr, sublists, sub_counts, comm_sz);

	//Gather partial list counts at root
	MPI_Gather(&count, 1, MPI_INT, recv_counts, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(my_rank == 0) {
		displacements[0] = 0;
		for(i = 1; i < comm_sz; ++i) {
			displacements[i] = displacements[i-1] + recv_counts[i-1];
		}
	}

	//Gather all partial lists at root
	MPI_Gatherv(my_arr, count, MPI_INT, arr, recv_counts, displacements, MPI_INT, 0, MPI_COMM_WORLD);

	if(my_rank == 0) {
		free(all_samples);
		free(recv_counts);
		free(displacements);
	}

	for(i = 0; i < comm_sz; ++i) {
		free(sublists[i]);
	}
	free(sublists);
	free(sub_counts);
	free(requests);
	free(pivots);
	free(samples);
	free(my_arr);
}

int partition(int arr[], int start, int end, int pivot) {
	int i;

	for(i = start; (i < end) && (arr[i] <= pivot); ++i);

	return i;
}

int merge(int arr[], int *sublists[], int list_counts[], int n_lists) {
	int *i_sublists = (int*)malloc(n_lists * sizeof(int)),
		*sub_values = (int*)malloc(n_lists * sizeof(int)),
		*value_valid = (int*)malloc(n_lists * sizeof(int));
	
	int i_arr = 0;
	
	int i;
	for(i = 0; i < n_lists; ++i) {
		i_sublists[i] = 0;
	}

	for(;;) {
		int sub_value_count = 0;
		for(i = 0; i < n_lists; ++i) {
			if(i_sublists[i] < list_counts[i]) {
				sub_values[i] = sublists[i][i_sublists[i]];
				value_valid[i] = 1;
				sub_value_count++;
			}
			else {
				value_valid[i] = 0;
			}
		}
		if(sub_value_count == 0) {
			break;
		}

		int min_list = min_index(sub_values, value_valid, n_lists);
		arr[i_arr++] = sublists[min_list][i_sublists[min_list]++];
	}

	free(i_sublists);
	free(sub_values);
	free(value_valid);

	return i_arr;
}

int min_index(int *values, int *mask, int n) {
	int min = INT_MAX, min_index = 0;

	int i;
	for(i = 0; i < n; ++i) {
		if(mask[i] && (values[i] < min)) {
			min = values[i];
			min_index = i;
		}
	}

	return min_index;
}
