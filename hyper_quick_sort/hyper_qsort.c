#include "hyper_qsort.h"
#include "serial_qsort.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

size_t hyper_qsort_rec(int arr[], int scratch[], int merge_scratch[], size_t size,
	size_t scratchSize, int blockStart, int blockEnd, int my_rank);

static size_t merge(int* in_result, size_t start, size_t stop, int* in_scratch,
	size_t scratchSize, int* merge_scratch); 
static size_t partition(int arr[], size_t size, int pivot);

static int hcube_level(int start, int end);
static int median(int arr[], int size);

static void print_array(int my_rank, int* arr, size_t size);

void hyper_qsort(int arr[], size_t size, int my_rank, int comm_sz) {
	int* scratch = (int*)malloc(size * sizeof(int));
	int* merge_scratch = (int*)malloc(size * sizeof(int));
	int* my_arr = (int*)malloc(size * sizeof(int));
	int count;
	int *recvCounts, *displacements;

	if(my_rank == 0) {
		printf("[Info] Starting list size: %d\n", size);

		recvCounts = (int*)malloc(comm_sz * sizeof(int));
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

	//Sort my array chunk using serial quicksort
	serial_qsort(my_arr, count);

	//Enter recursive hyper_qsort routine
	/*
	int level;
	for(level = max_level(comm_sz); level >= 0; --level) {
		count = hyper_qsort_rec(my_arr, scratch, merge_scratch, count, size, level, my_rank,
			comm_sz);
	}
	*/
	count = hyper_qsort_rec(my_arr, scratch, merge_scratch, count, size, 0, comm_sz, my_rank);
	
	MPI_Gather(&count, 1, MPI_INT, recvCounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(my_rank == 0) {
		int i;
		displacements[0] = 0;
		for(i = 1; i < comm_sz; ++i) {
			displacements[i] = displacements[i-1] + recvCounts[i-1];
		}
	}
	MPI_Gatherv(my_arr, count, MPI_INT, arr, recvCounts, displacements, MPI_INT,
		0, MPI_COMM_WORLD);
	
	free(my_arr);
	free(scratch);
}


size_t hyper_qsort_rec(int arr[], int scratch[], int merge_scratch[], size_t size,
	size_t scratchSize, int blockStart, int blockEnd, int my_rank) {
	
	if((blockEnd - blockStart) < 2) {
		//End of recursion
		printf("[%d] Block (%d, %d): End of recursion (my values in range [%d, %d])\n", my_rank,
			blockStart, blockEnd, arr[0], arr[size-1]);
		return size;
	}

	printf("[%d] Initial list size: %d\n", my_rank, size);
	
	MPI_Status status;
	int split = blockStart + (blockEnd - blockStart) / 2 + ((blockEnd - blockStart) % 2),
		block_rank = my_rank - blockStart, lowerSubBlockSize = (split - blockStart),
		upperSubBlockSize = (blockEnd - split), subBlockStart, subBlockEnd;
	int pivot;

	if(my_rank == blockStart) {
		//Root always provides the pivot
		pivot = median(arr, size);
		
		//printf("Level=%d, pivot=%d\n", level, pivot);
		printf("[%d] Sub-block sizes: %d, %d\n", my_rank, lowerSubBlockSize, upperSubBlockSize);
		printf("[%d] Sending pivot to group (%d, %d)\n", my_rank, blockStart, blockEnd);
		int i;
		for(i = blockStart+1; i < blockEnd; ++i) {
			MPI_Send(&pivot, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		}
	}
	else {
		printf("[%d] Waiting for pivot in group (%d, %d)\n", my_rank, blockStart, blockEnd);
		MPI_Recv(&pivot, 1, MPI_INT, blockStart, 0, MPI_COMM_WORLD, &status);
	}

	printf("[%d] Pivot exchange complete\n", my_rank);

	size_t i_pivot = partition(arr, size, pivot);

	if(my_rank < split) {
		subBlockStart = blockStart;
		subBlockEnd = split;

		//Calculating the upper neighbor is more complicated if block sizes can be non powers of two
		//We mod the optimal (power of 2) neighbor with the actual upper sub-block size to efficiently
		//Split the work with the actual number of available upper block processes
		int neighbor = (block_rank % upperSubBlockSize) + split;
		int recv_count;

		printf("[%d] Single neighbor is %d\n", my_rank, neighbor);
		printf("[%d] Sending %d values to neighbor\n", my_rank, size - i_pivot);

		//Send upper list to neighbor
		MPI_Send(arr + i_pivot, size - (i_pivot), MPI_INT, neighbor, 0,
			MPI_COMM_WORLD);

		//Receive neighbor's lower list
		MPI_Recv(scratch, scratchSize, MPI_INT, neighbor, 0, MPI_COMM_WORLD,
			&status);
		MPI_Get_count(&status, MPI_INT, &recv_count);

		printf("[%d] Received %d values from neighbor\n", my_rank, recv_count);

		//Merge lists into sorted intermediate result
		size = merge(arr, 0, i_pivot, scratch, recv_count, merge_scratch);
		if(!validate(arr, size)) {
			printf("[%d] Validation error after merge\n", my_rank);
		}
	}
	else {
		subBlockStart = split;
		subBlockEnd = blockEnd;
		int subBlockRank = my_rank - subBlockStart;

		//Because block sizes can be non powers of two, each upper sub-block process may have
		//more than one neighbor
		int neighbor_count = (lowerSubBlockSize / upperSubBlockSize) +
			(((lowerSubBlockSize % upperSubBlockSize) > subBlockRank) ? 1 : 0);
		int sendSize = i_pivot;
		int scratchEnd = 0;

		printf("[%d] Neighbor count: %d (%d, %d)\n", my_rank, neighbor_count, lowerSubBlockSize,
			upperSubBlockSize);

		int* partialList = (int*)malloc(scratchSize * sizeof(int));

		int i;
		for(i = 0; i < neighbor_count; ++i) {
			int neighbor = blockStart + subBlockRank + i*upperSubBlockSize;

			printf("[%d] Swapping lists with neighbor %d\n", my_rank, neighbor);

			//Receive this neighbor's upper list
			int recv_count;
			MPI_Recv((scratchEnd > 0) ? partialList : scratch, scratchSize, MPI_INT, neighbor, 0,
				MPI_COMM_WORLD, &status);
			MPI_Get_count(&status, MPI_INT, &recv_count);
			printf("[%d] Received %d values from neighbor\n", my_rank, recv_count);
			
			if(scratchEnd > 0) {
				//Merge this upper list with already received upper lists
				scratchEnd = merge(scratch, 0, scratchEnd, partialList, recv_count, merge_scratch);
			}
			else {
				scratchEnd = recv_count;
			}

			//Send part of lower list to this neighbor
			int sendStart = i*sendSize/neighbor_count, sendEnd = (i+1)*sendSize/neighbor_count;
			MPI_Send(arr + sendStart, sendEnd - sendStart, MPI_INT, neighbor, 0, MPI_COMM_WORLD);
			printf("[%d] Sending %d values to neighbor\n", my_rank, sendEnd - sendStart);
		}

		//Merge all received lists with my current list
		size = merge(arr, i_pivot, size, scratch, scratchEnd, merge_scratch);
		if(!validate(arr, size)) {
			printf("[%d] Validation error after merge\n", my_rank);
		}

		free(partialList);
	}

	printf("[%d] List size: %d\n", my_rank, size);
	
	printf("[%d] Splitting into block (%d, %d)\n", my_rank, subBlockStart, subBlockEnd);
	size = hyper_qsort_rec(arr, scratch, merge_scratch, size, scratchSize, subBlockStart,
		subBlockEnd, my_rank);

	return size;
}

size_t merge(int* in_result, size_t start, size_t stop, int* in_scratch,
	size_t scratchSize, int* merge_scratch) {

	if(!validate(in_result, stop-start)) {
		printf("[Error] merge: in_result not sorted\n");
	}
	if(!validate(in_scratch, scratchSize)) {
		printf("[Error] merge: scratch not sorted\n");
	}

	size_t i_out, i_result = start, i_scratch = 0;
	for(i_out = 0; (i_result < stop) && (i_scratch < scratchSize); ++i_out) {
		if(in_result[i_result] <= in_scratch[i_scratch]) {
			//Move value from in_result list to output list
			merge_scratch[i_out] = in_result[i_result];
			++i_result;
		}
		else {
			//Move value from scratch list to output list
			merge_scratch[i_out] = in_scratch[i_scratch];
			i_scratch++;
		}
	}
	if(!validate(merge_scratch, i_out)) {
		printf("[Error] merge: output before final pass not sorted\n");
	}
	//Copy any leftover values
	for(; i_result < stop; ++i_result, ++i_out) {
		merge_scratch[i_out] = in_result[i_result];
	}
	if(!validate(merge_scratch, i_out)) {
		printf("[Error] merge: output before scratch pass not sorted\n");
	}
	for(; i_scratch < scratchSize; ++i_scratch, ++i_out) {
		merge_scratch[i_out] = in_scratch[i_scratch];
	}

	if(!validate(merge_scratch, i_out)) {
		printf("[Error] merge: output not sorted\n");
	}

	memcpy(in_result, merge_scratch, i_out * sizeof(int));
	if(!validate(in_result, i_out)) {
		printf("[Error] merge: output not sorted\n");
	}

	return i_out;
}

size_t partition(int arr[], size_t size, int pivot) {
	size_t i;
	for(i = 0; (i < size) && (arr[i] <= pivot); ++i);

	return i;
}

int hcube_level(int start, int end) {
	int level = 0;
	int size = end - start - 1;

	for(; size > 0; size >>= 1) {
		++level;
	}

	return level - 1;
}

int median(int arr[], int size) {
	if(size & 0x01) {
		//Odd number size
		return (arr[size/2] + arr[size/2 + 1])/2;
	}
	else {
		//Even number size
		return arr[size/2];
	}
}

void print_array(int my_rank, int* arr, size_t size) {
	int i;
	for(i = 0; i < size; ++i) {
		printf("[%d] %d\n", my_rank, arr[i]);
	}
	printf("\n");
}
