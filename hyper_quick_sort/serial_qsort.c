#include "serial_qsort.h"

#include <stdio.h>

void serial_qsort_rec(int* arr, size_t start, size_t stop);
static size_t partition(int* arr, size_t start, size_t stop);

void serial_qsort(int* arr, size_t size) {
  serial_qsort_rec(arr, 0, size-1);
}

void serial_qsort_rec(int* arr, size_t start, size_t stop) {
  if(start < stop) {
    size_t p = partition(arr, start, stop);
    if(p > 0) {
      serial_qsort_rec(arr, start, p-1);
    }
    if(p < stop) {
      serial_qsort_rec(arr, p+1, stop);
    }
  }
}

size_t partition(int* arr, size_t start, size_t stop) {
  int pivot = arr[stop];
  
  ssize_t i = start-1, j;
  for(j = start; j < stop; ++j) {
    if(arr[j] <= pivot) {
      ++i;
      swap(&arr[i], &arr[j]);
    }
  }
  swap(&arr[i+1], &arr[stop]);
  
  return i+1;
}

int validate(int* arr, size_t size) {
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

void swap(int* a, int* b) {
	int t = *a;
	*a = *b;
	*b = t;
}
