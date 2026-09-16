#ifndef RAPTOR_UTILITY_H
#define RAPTOR_UTILITY_H

void align_xor(unsigned char *dst, unsigned char *src, unsigned int len);


unsigned int align_xor_pop(uint64_t *dst, uint64_t *src, unsigned int cnt);

unsigned int align_test_include(uint64_t *mask, uint64_t *data, unsigned int cnt);


int compare_uint32(const void *elem1, const void *elem2);

void selection_sort(unsigned int *list, unsigned int max);

void bubble_sort(unsigned int *list, unsigned int max);


void print_matrix(
	unsigned char *matrix,
	unsigned int num_col,	// number of columns = size of row
	unsigned int num_row);	// number of rows = size of column

#endif
