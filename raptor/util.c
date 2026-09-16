/*
MIT License
Copyright (c) 2026 Yutaka Sawada
*/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <memory.h>

/*
// To use builtin popcount for MSVC
#ifdef _MSC_VER
#include <intrin.h>
#endif
*/


// Both dst and src must be aligned.
void align_xor(unsigned char *dst, unsigned char *src, unsigned int len)
{
	len /= 8;	// convert from bytes to number of uint64_t.
	for (unsigned int i = 0; i < len; i++)
		((uint64_t *)dst)[i] ^= ((uint64_t *)src)[i];
}

// XOR and popcount
#define M1 0x5555555555555555
#define M2 0x3333333333333333
#define M4 0x0f0f0f0f0f0f0f0f
unsigned int align_xor_pop(
	uint64_t *dst, uint64_t *src,
	unsigned int cnt)	// number of uint64_t
{
	uint32_t i, pops = 0;
	uint64_t x;

	for (i = 0; i < cnt; i++){
		x = dst[i] ^ src[i];
		dst[i] = x;

		// popcount64b from https://en.wikipedia.org/wiki/Hamming_weight
		x -= (x >> 1) & M1;
		x = (x & M2) + ((x >> 2) & M2);
		x = (x + (x >> 4)) & M4;
		x += x >>  8;
		x += x >> 16;
		x += x >> 32;
		pops += x & 0x7f;

		// speed difference is small
		// pops += (uint32_t)__popcnt64(x);
	}

	return pops;
}

// When mask bits includes all bits in data, it returns -1.
// Or else, it returns position of the first different bit.
unsigned int align_test_include(
	uint64_t *mask,
	uint64_t *data,
	unsigned int cnt)	// number of uint64_t
{
	uint32_t i;
	uint64_t mask2, data2;

	for (i = 0; i < cnt; i++){
		mask2 = mask[i];
		data2 = data[i];
		mask2 = mask2 & data2;
		if (mask2 != data2){
			uint32_t pops;
			uint64_t x = data2 ^ mask2;
			x = (x & (-(int64_t)x)) - 1;

			// popcount64b
			x -= (x >> 1) & M1;
			x = (x & M2) + ((x >> 2) & M2);
			x = (x + (x >> 4)) & M4;
			x += x >>  8;
			x += x >> 16;
			x += x >> 32;
			pops = x & 0x7f;

			return i * 64 + pops;
		}
	}

	return 0xffffffff;
}


// compare function for C runtime library
int compare_uint32(const void *elem1, const void *elem2)
{
	unsigned int val1 = *((unsigned int *)elem1);
	unsigned int val2 = *((unsigned int *)elem2);

	if (val1 < val2)
		return -1;
	if (val1 > val2)
		return 1;
	return 0;
}

void selection_sort(unsigned int *list, unsigned int max)
{
	unsigned int i, j, min, i_min;

	for (i = 0; i < max - 1; i++){
		min = list[i];
		i_min = i;
		for (j = i + 1; j < max; j++){
			if (min > list[j]){
				min = list[j];
				i_min = j;
			}
		}
		if (i_min != i){
			list[i_min] = list[i];
			list[i] = min;
		}
	}
}

void bubble_sort(unsigned int *list, unsigned int max)
{
	unsigned int i, tmp, range;

	for (range = max; range > 1; range--){
		for (i = 0; i < range - 1; i++){
			if (list[i] > list[i + 1]){
				tmp = list[i];
				list[i] = list[i + 1];
				list[i + 1] = tmp;
			}
		}
	}
}


void print_matrix(
	unsigned char *matrix,
	unsigned int num_col,	// number of columns = size of row
	unsigned int num_row)	// number of rows = size of column
{
	unsigned int col_id, row_id;

	printf("matrix %u * %u\n", num_col, num_row);
	for (row_id = 0; row_id < num_row; row_id++){
		for (col_id = 0; col_id < num_col; col_id++){
			printf(" %d", matrix[num_col * row_id + col_id]);
		}
		printf("\n");
	}
}



