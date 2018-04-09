#include "sort_util.h"

int validate(int arr[], size_t size) {
	if(size < 2) {
		return 1;
	}

	size_t i;
	for(i = 0; i < (size-1); ++i) {
		if(arr[i] > arr[i+1]) {
			return 0;
		}
	}

	return 1;
}
