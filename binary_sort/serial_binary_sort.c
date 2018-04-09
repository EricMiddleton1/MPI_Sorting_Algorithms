#include "serial_binary_sort.h"

#include <string.h>
#include <stdio.h>

static int binary_search(int arr[], size_t start, size_t stop, int value);

void serial_binary_sort(int arr[], size_t size) {
	int i;
	for(i = 1; i < size; ++i) {
		int value = arr[i];
		int insert_loc = binary_search(arr, 0, i, value);

		if(insert_loc < i) {
			memmove(arr + insert_loc + 1, arr + insert_loc, (i - insert_loc) * sizeof(int));
		}
		arr[insert_loc] = value;
	}
}

int binary_search(int arr[], size_t start, size_t stop, int value) {
	size_t middle = start + (stop - start)/2;

	if((stop - start) <= 1) {
		//Base case
		return (value > arr[start]) ? (start + 1) : start;
	}

	if(value < arr[middle]) {
		//Recursively search through lower half
		return binary_search(arr, start, middle, value);
	}
	else {
		//Recursively search through upper half
		return binary_search(arr, middle, stop, value);
	}
}
