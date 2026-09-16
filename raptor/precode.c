/*
MIT License
Copyright (c) 2026 Yutaka Sawada
*/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "raptor.h"
#include "precode.h"
#include "bitmask.h"
#include "gene.h"
#include "util.h"


// 5.4.2.3. Pre-Coding Relationships

// function from nanorq
static int is_prime(int n) {
	if (n <= 1)
		return 0;
	if (n <= 3)
		return 1;

	if (n % 2 == 0 || n % 3 == 0)
		return 0;

	for (int i = 5; i * i <= n; i = i + 6)
		if (n % i == 0 || n % (i + 2) == 0)
			return 0;

	return 1;
}

// the number of ways j objects can be chosen from among i objects without repetition.
static unsigned int choose(unsigned int i, unsigned int j){
	unsigned int ci, cj; // combination
	unsigned int count;

	ci = i;
	for (count = 1; count < j; count++)
		ci *= i - count;
	cj = j;
	for (count = 1; count < j; count++)
		cj *= j - count;

	//printf("i = %u, j = %u, ci = %u, cj = %u, ci/cj = %u\n", i, j, ci, cj, ci / cj);
	return ci / cj;
}

void init_precode(raptor_coder *rc)
{
	// X be the smallest positive integer such that X*(X-1) >= 2*K.
	unsigned int K = rc->K;
	unsigned int X = 4;
	// Because minimum 2*K = 4 * 2 = 8, minimum X = 4.
	// Because maximum 2*K = 8192 * 2 = 16384, maximum X = 129.
	while (X * (X - 1) < 2 * K)
		X++;
	//printf("X = %u\n", X);

	// S be the smallest prime integer such that S >= ceil(0.01*K) + X
	unsigned int S = (K + 99) / 100 + X;
	while (is_prime(S) == 0)
    	S++;
	rc->S = S;
	// Minimum S = 1 + 4 = 5
	// Maximun S = 82 + 129 = 211

	// H be the smallest integer such that choose(H,ceil(H/2)) >= K + S
	unsigned int H = 5;
	// Minimum K + S = 4 + 5 = 9, H = 5
	// Maximun K + S = 8192 + 211, H = 16
	while (choose(H, (H + 1) / 2) < K + S)
		H++;
	rc->H = H;

	// L = K+S+H
	// Maximun L = 8419, Ld = 8419
	unsigned int L = K + S + H;
	rc->L = L;

	// L' be the smallest prime that is greater than or equal to L
	unsigned int Ld = L;
	while (is_prime(Ld) == 0)
    	Ld++;
	rc->Ld = Ld;

	printf("K = %u, S = %u, H = %u, L = %u, Ld = %u\n", K, S, H, L, Ld);
}


/*
To make these functions, I refer to "gofountain".

gofountain
Copyright 2014 Google Inc.
Apache License Version 2.0
*/

// bitsSet returns how many bits in x are set.
// This algorithm basically uses shifts and ANDs to sum up the bits in
// a tree fashion.
static unsigned int bitsSet(uint64_t x){
	x -= (x >> 1) & 0x5555555555555555;
	x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333);
	x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0f;
	return (unsigned int)((x * 0x0101010101010101) >> 56);
}

// grayCode calculates the gray code representation of the input argument
// The Gray code is a binary representation in which successive values differ
// by exactly one bit. See http://en.wikipedia.org/wiki/Gray_code
static uint64_t grayCode(uint64_t x){
	return (x >> 1) ^ x;
}

// buildGraySequence returns a sequence (in ascending order) of "length" Gray numbers,
// all of which have exactly "b" bits set.
void buildGraySequence(unsigned int *seq, unsigned int length, unsigned int b){
	unsigned int i = 0;
	uint64_t x = 0;
	uint64_t g;
	while (i < length){
		g = grayCode(x);
		if (bitsSet(g) == b){
			seq[i] = (unsigned int)g;
			i++;
		}
		x++;
	}
}

// returns number of 1s in matrix, 0 = fail
unsigned int init_precode_matrix(
	raptor_coder *rc,
	unsigned char *A)		// zero fill at first
{
	unsigned int i, a, b, h, j;
	unsigned int K = rc->K;
	unsigned int S = rc->S;
	unsigned int H = rc->H;
	unsigned int Hd = (H + 1) / 2;
	unsigned int L = rc->L;		// number of columns in matrix A
	unsigned int *m;
	unsigned int count;

	// G_LDPC be the S x K generator matrix of the LDPC symbols.
	for (i = 0; i < K; i++){
		a = 1 + ((i / S) % (S - 1));
		b = i % S;
		A[L * b + i] = 1;
		b = (b + a) % S;
		A[L * b + i] = 1;
		b = (b + a) % S;
		A[L * b + i] = 1;
	}
	//printf("G_LDPC S x K = %u x %u\n", S, K);
	count = K * 3;	// column weight 3

	// G_Half be the H x (K+S) generator matrix of the Half symbols
	m = malloc(sizeof(unsigned int) * (K + S));
	if (m == NULL)
		return 0;
	buildGraySequence(m, K + S, Hd);
	for (h = 0; h < H; h++){
		for (j = 0; j < K + S; j++){
			if ((m[j] >> h) & 1){
				A[L * (S + h) + j] = 1;
				count++;
			}
		}
	}
	free(m);
	//printf("G_Half H x (K+S) = %u x %u\n", H, K + S);

	// I_S be the S x S identity matrix
	for (i = 0; i < S; i++)
		A[L * i + K + i] = 1;
	count += S;

	// I_H be the H x H identity matrix
	for (i = 0; i < H; i++)
		A[L * (S + i) + K + S + i] = 1;
	count += H;

	return count;
}

// returns number of 1s in matrix, 0 = fail
int init_precode_bit_matrix(
	raptor_coder *rc,
	uint64_t *A,			// zero fill at first
	unsigned int int_col)	// number of integers for columns in matrix A
{
	unsigned int i, a, b, h, j;
	unsigned int K = rc->K;
	unsigned int S = rc->S;
	unsigned int H = rc->H;
	unsigned int Hd = (H + 1) / 2;
	unsigned int L = rc->L;
	unsigned int *m;
	unsigned int count;

	// G_LDPC be the S x K generator matrix of the LDPC symbols.
	for (i = 0; i < K; i++){
		a = 1 + ((i / S) % (S - 1));
		b = i % S;
		bitmask_set(A + int_col * b, i);
		b = (b + a) % S;
		bitmask_set(A + int_col * b, i);
		b = (b + a) % S;
		bitmask_set(A + int_col * b, i);
	}
	//printf("G_LDPC S x K = %u x %u\n", S, K);
	count = K * 3;	// column weight 3

	// G_Half be the H x (K+S) generator matrix of the Half symbols
	m = malloc(sizeof(unsigned int) * (K + S));
	if (m == NULL)
		return 0;
	buildGraySequence(m, K + S, Hd);
	for (h = 0; h < H; h++){
		for (j = 0; j < K + S; j++){
			if ((m[j] >> h) & 1){
				bitmask_set(A + int_col * (S + h), j);
				count++;
			}
		}
	}
	free(m);
	//printf("G_Half H x (K+S) = %u x %u\n", H, K + S);

	// I_S be the S x S identity matrix
	for (i = 0; i < S; i++)
		bitmask_set(A + int_col * i, K + i);
	count += S;

	// I_H be the H x H identity matrix
	for (i = 0; i < H; i++)
		bitmask_set(A + int_col * (S + i), K + S + i);
	count += H;

	return count;
}

// returns 0 = success, others = fail
int make_generator_matrix(raptor_coder *rc)
{
	unsigned int i, a, b, d, j;
	unsigned int K = rc->K;
	unsigned int S = rc->S;
	unsigned int H = rc->H;
	unsigned int L = rc->L;
	unsigned int Ld = rc->Ld;
	unsigned int int_col = (L + IDXBITS - 1) / IDXBITS;	// number of integers for columns
	unsigned short *order_list, *degree_list;
	uint64_t *A, *row_p;

	// allocate two lists at once
	order_list = malloc(sizeof(unsigned short) * L * 2);
	if (order_list == NULL)
		return 1;
	degree_list = order_list + L;

	// allocate matrix of L * L by using integers
	A = calloc(int_col * L, sizeof(uint64_t));
	if (A == NULL){
		free(order_list);
		return 2;
	}
	//printf("bit matrix's %u columns -> %u integers, %g KB\n", L, int_col, (double)(int_col * L) / 125.0);

	// make constraint matrix of the precode
	if (init_precode_bit_matrix(rc, A, int_col) == 0){
		free(order_list);
		free(A);
		return 3;
	}
	//print_bit_matrix(A, L, S + H);

	// set degree of check symbols
	row_p = A;
	for (i = 0; i < S + H; i++){
		degree_list[i] = bitmask_pop(row_p, int_col);
		row_p += int_col;
	}

	// G_LT be the KxL generator matrix of the encoding symbols generated by the LT Encoder.
	row_p = A + int_col * (S + H);
	for (i = 0; i < K; i++){	// ESI
		// 5.4.4.3.  LT Encoding Symbol Generator
		Trip(K, Ld, i, &d, &a, &b);
		//printf("d = %u, a = %u, b = %u\n", d, a, b);
		if (d > L)
			d = L;
		while (b >= L)
			b = (b + a) % Ld;
		bitmask_set(row_p, b);
		for (j = 1; j < d; j++){
			b = (b + a) % Ld;
			while (b >= L)
				b = (b + a) % Ld;
			bitmask_set(row_p, b);
		}

		degree_list[S + H + i] = d;	// set degree of encoding symbol
		row_p += int_col;
	}
	//print_bit_matrix(A, L, L);

	// solve equation
	unsigned int col_id, row_id, pivot_id, next_id;
	unsigned int row_degree, pivot_degree;
	unsigned int Tal = rc->Tal;
	unsigned char *Cd = rc->tmp_buf;	// array of source symbols (S+H ~ L-1)
	unsigned char *C_p;					// pointer of intermediate symbols
	uint64_t *pivot_p;

	// find a row with the smallest degree
	pivot_degree = 0xffff;
	for (i = 0; i < L; i++){
		row_degree = degree_list[i];
		//printf("degree_list[%u] = %u\n", i, row_degree);
		if (pivot_degree > row_degree){
			pivot_degree = row_degree;
			pivot_id = i;
		}
	}

	//a = 0;
	while (pivot_degree < 0xffff){
		col_id = bitmask_ntz(A + int_col * pivot_id, int_col);
		/*
		if (col_id >= L){	// Because encoder should not fail, no need this error check.
			printf("Cannot solve equation, pivot_id = %u, col_id = %u\n", pivot_id, col_id);
			free(order_list);
			free(A);
			return 4;
		}
		*/
		//printf("pivot_id = %u, pivot_degree = %u, col_id = %u\n", pivot_id, pivot_degree, col_id);
		pivot_p = A + int_col * pivot_id;
		order_list[pivot_id] = col_id;	// save order of rows
		degree_list[pivot_id] = 0xffff;	// no need degree of pivot row

		// XOR pivot_row to other rows
		pivot_degree = 0xffff;
		C_p = Cd + (size_t)Tal * pivot_id;
		row_p = A;
		for (row_id = 0; row_id < L; row_id++){
			if (row_id == pivot_id){
				row_p += int_col;
				continue;
			}

			if (bitmask_check(row_p, col_id) != 0){	// XOR rows
				align_xor(Cd + (size_t)Tal * row_id, C_p, Tal);
//				for (i = 0; i < int_col; i++)
//					row_p[i] ^= pivot_p[i];
				if (degree_list[row_id] != 0xffff){	// re-calculate degree
					//row_degree = bitmask_pop(row_p, int_col);
					row_degree = align_xor_pop(row_p, pivot_p, int_col);
					/*
					if (row_degree == 0){	// If row_degree is zero, this function would fail.
						printf("Cannot solve equation, pivot_id = %u\n", pivot_id);
						free(order_list);
						free(A);
						return 5;
					}
					*/
					degree_list[row_id] = row_degree;
					//printf("degree_list[%u] = %u\n", row_id, row_degree);
					if (pivot_degree > row_degree){
						pivot_degree = row_degree;
						next_id = row_id;
					}
				} else {
					for (i = 0; i < int_col; i++)
						row_p[i] ^= pivot_p[i];
				}

			} else if (degree_list[row_id] != 0xffff){	// find a row with the smallest degree
				row_degree = degree_list[row_id];
				if (pivot_degree > row_degree){
					pivot_degree = row_degree;
					next_id = row_id;
				}
			}
			row_p += int_col;
		}

/*
		if (a >= 3)
			print_bit_matrix(A, L, L);
		if (a >= 5)
			break;
		a++;
*/
		pivot_id = next_id;	// set next pivot
	}
	//print_bit_matrix(A, L, L);

	// Because encoder doesn't fail, no need matrix anymore.
	free(A);

	// restore intermediate symbols
	C_p = rc->C;
	for (row_id = 0; row_id < L; row_id++){
		col_id = order_list[row_id];	// load order of rows
		memcpy(C_p + (size_t)Tal * col_id, Cd, Tal);

		Cd += Tal;	// goto next symbol
	}
	free(order_list);

	return 0;
}

/*
refer to a papar by Francisco Lazaro and Gerhard Bauch;
Inactivation Decoding of LT and Raptor Codes: Analysis and Code Design

returns 0 = success, others = fail
*/
int make_inactive_matrix(raptor_coder *rc)
{
	unsigned int i, a, b, d, j;
	unsigned int K = rc->K;
	unsigned int S = rc->S;
	unsigned int H = rc->H;
	unsigned int L = rc->L;
	unsigned int Ld = rc->Ld;
	unsigned int int_col = (L + IDXBITS - 1) / IDXBITS;	// number of integers for columns
	unsigned short *order_list, *degree_list;
	uint64_t *A, *row_p;
	uint64_t *act_mask;

	// allocate two lists at once
	order_list = malloc(sizeof(unsigned short) * L * 2);
	if (order_list == NULL)
		return 1;
	degree_list = order_list + L;

	// setup list of intermediate symbols
	act_mask = (uint64_t *)calloc(int_col, sizeof(uint64_t));
	if (act_mask == NULL){
		free(order_list);
		return 2;
	}

	// allocate matrix of L * L by using integers
	A = calloc(int_col * L, sizeof(uint64_t));
	if (A == NULL){
		free(order_list);
		free(act_mask);
		return 3;
	}
	//printf("bit matrix's %u columns -> %u integers, %g KB\n", L, int_col, (double)(int_col * L) / 125.0);

	// make constraint matrix of the precode
	if (init_precode_bit_matrix(rc, A, int_col) == 0){
		free(order_list);
		free(act_mask);
		free(A);
		return 4;
	}
	//print_bit_matrix(A, L, S + H);

	// set degree of check symbols
	row_p = A;
	for (i = 0; i < S + H; i++){
		degree_list[i] = bitmask_pop(row_p, int_col);
		row_p += int_col;
	}

	// G_LT be the KxL generator matrix of the encoding symbols generated by the LT Encoder.
	row_p = A + int_col * (S + H);
	for (i = 0; i < K; i++){	// ESI
		// 5.4.4.3.  LT Encoding Symbol Generator
		Trip(K, Ld, i, &d, &a, &b);
		//printf("d = %u, a = %u, b = %u\n", d, a, b);
		if (d > L)
			d = L;
		while (b >= L)
			b = (b + a) % Ld;
		bitmask_set(row_p, b);
		for (j = 1; j < d; j++){
			b = (b + a) % Ld;
			while (b >= L)
				b = (b + a) % Ld;
			bitmask_set(row_p, b);
		}

		degree_list[S + H + i] = d;	// set degree of encoding symbol
		row_p += int_col;
	}
	//print_bit_matrix(A, L, L);

	// solve equation
	unsigned int col_id, row_id, pivot_id;
	unsigned int next_pivot, next_col;
	unsigned int row_degree, pivot_degree;
	unsigned int Tal = rc->Tal;
	unsigned char *Cd = rc->tmp_buf;	// array of source symbols (S+H ~ L-1)
	unsigned char *C_p;					// pointer of intermediate symbols
	uint64_t *pivot_p;

	// find a row with the smallest degree
	pivot_degree = 0xfffd;
	for (i = 0; i < L; i++){
		row_degree = degree_list[i];
		//printf("degree_list[%u] = %u\n", i, row_degree);
		if (pivot_degree > row_degree){
			pivot_degree = row_degree;
			pivot_id = i;
		}
	}
	col_id = bitmask_ntz(A + int_col * pivot_id, int_col);

	// (a) Triangulation process
	// (b) Zero matrix procedure
	//a = 0;
	while (pivot_degree < 0xfffd){
		/*
		if (col_id >= L){	// Because encoder should not fail, no need this error check.
			printf("Cannot solve equation, pivot_id = %u, col_id = %u\n", pivot_id, col_id);
			free(order_list);
			free(act_mask);
			free(A);
			return 5;
		}
		*/
		//printf("pivot_id = %u, pivot_degree = %u, col_id = %u\n", pivot_id, pivot_degree, col_id);
		pivot_p = A + int_col * pivot_id;
		order_list[pivot_id] = col_id;	// save order of rows
		degree_list[pivot_id] = 0xffff;	// no need degree of pivot row

		if (pivot_degree >= 2){	// inactivate following symbols
			for (i = col_id / IDXBITS; i < int_col; i++)
				act_mask[i] |= pivot_p[i];
			//print_bit_matrix(act_mask, L, 1);
		}

		// XOR pivot_row to other rows
		pivot_degree = 0xfffd;
		C_p = Cd + (size_t)Tal * pivot_id;
		row_p = A;
		for (row_id = 0; row_id < L; row_id++){
			if (degree_list[row_id] >= 0xfffd){	// skip old pivot rows and inactive rows
				//printf("skip row %u, degree_list = %u\n", row_id, degree_list[row_id]);
				row_p += int_col;
				continue;
			}

			if (bitmask_check(row_p, col_id) != 0){	// XOR rows
				align_xor(Cd + (size_t)Tal * row_id, C_p, Tal);
//				for (i = 0; i < int_col; i++)
//					row_p[i] ^= pivot_p[i];
				// re-calculate degree
//				row_degree = bitmask_pop(row_p, int_col);
				row_degree = align_xor_pop(row_p, pivot_p, int_col);
				/*
				if (row_degree == 0){	// If row_degree is zero, this function would fail.
					printf("Cannot solve equation, pivot_id = %u\n", pivot_id);
					free(order_list);
					free(act_mask);
					free(A);
					return 6;
				}
				*/
				if (pivot_degree > row_degree){
					// test inactive row
					j = align_test_include(act_mask, row_p, int_col);
					if (j == 0xffffffff){	// all symbols are inactive
						row_degree = 0xfffd;
						//printf("inactivate XORed row %u\n", row_id);
					} else {
						pivot_degree = row_degree;
						next_pivot = row_id;
						next_col = j;
					}
				}
				degree_list[row_id] = row_degree;
				//printf("degree_list[%u] = %u\n", row_id, row_degree);

			} else {	// find a row with the smallest degree
				row_degree = degree_list[row_id];
				if (pivot_degree > row_degree){
					// test inactive row
					j = align_test_include(act_mask, row_p, int_col);
					if (j == 0xffffffff){	// all symbols are inactive
						degree_list[row_id] = 0xfffd;
						//printf("inactivate row %u\n", row_id);
					} else {
						pivot_degree = row_degree;
						next_pivot = row_id;
						next_col = j;
					}
				}
			}
			row_p += int_col;
		}

/*
		printf("next_pivot = %u, next_col = %u, pivot_degree = %u\n", next_pivot, next_col, pivot_degree);
		if (a >= 2)
			print_bit_matrix(A, L, L);
		if (a >= 25)
			break;
		a++;
*/
		pivot_id = next_pivot;	// set next pivot
		col_id = next_col;
	}
	free(act_mask);	// No need to check active symbols anymore
	//print_bit_matrix(A, L, L);

	// find a row with the smallest degree
	//b = 0;	// count number of inactive rows
	pivot_degree = 0xfffd;
	row_p = A;
	for (i = 0; i < L; i++){
		row_degree = degree_list[i];
		if (row_degree == 0xfffd){	// re-calculate degree
			row_degree = bitmask_pop(row_p, int_col);
			degree_list[i] = row_degree;
			//b++;
		}
		//printf("degree_list[%u] = %u, order_list = %u\n", i, row_degree, order_list[i]);
		if (pivot_degree > row_degree){
			pivot_degree = row_degree;
			pivot_id = i;
		}
		row_p += int_col;
	}
	//printf("number of resolvable rows = %u, inactive rows = %u\n", L - b, b);

	// (c) Gaussian elimination
	while (pivot_degree < 0xfffd){
		//b--;
		col_id = bitmask_ntz(A + int_col * pivot_id, int_col);
		/*
		if (col_id >= L){	// Because encoder should not fail, no need this error check.
			printf("Cannot solve equation, pivot_id = %u, col_id = %u\n", pivot_id, col_id);
			free(order_list);
			free(A);
			return 7;
		}
		*/
		//printf("pivot_id = %u, pivot_degree = %u, col_id = %u\n", pivot_id, pivot_degree, col_id);
		pivot_p = A + int_col * pivot_id;
		order_list[pivot_id] = col_id;	// save order of rows
		degree_list[pivot_id] = 0xfffd;	// pivot row of inactive symbols

		// XOR pivot_row to other rows
		pivot_degree = 0xfffd;
		C_p = Cd + (size_t)Tal * pivot_id;
		row_p = A;
		for (row_id = 0; row_id < L; row_id++){
			if ((row_id == pivot_id) || (degree_list[row_id] == 0xffff)){	// skip old pivot rows
				row_p += int_col;
				continue;
			}

			if (bitmask_check(row_p, col_id) != 0){	// XOR rows
				align_xor(Cd + (size_t)Tal * row_id, C_p, Tal);
//				for (i = 0; i < int_col; i++)
//					row_p[i] ^= pivot_p[i];
				if (degree_list[row_id] != 0xfffd){	// re-calculate degree
					//row_degree = bitmask_pop(row_p, int_col);
					row_degree = align_xor_pop(row_p, pivot_p, int_col);
					/*
					if (row_degree == 0){	// If row_degree is zero, this function would fail.
						//printf("Cannot solve equation, pivot_id = %u\n", pivot_id);
						free(order_list);
						free(A);
						return 8;
					}
					*/
					degree_list[row_id] = row_degree;
					//printf("degree_list[%u] = %u\n", row_id, row_degree);
					if (pivot_degree > row_degree){
						pivot_degree = row_degree;
						next_pivot = row_id;
					}
				} else {
					for (i = 0; i < int_col; i++)
						row_p[i] ^= pivot_p[i];
				}

			} else if (degree_list[row_id] != 0xfffd){	// find a row with the smallest degree
				row_degree = degree_list[row_id];
				if (pivot_degree > row_degree){
					pivot_degree = row_degree;
					next_pivot = row_id;
				}
			}
			row_p += int_col;
		}

		//printf("next_pivot = %u, pivot_degree = %u\n", next_pivot, pivot_degree);
		pivot_id = next_pivot;	// set next pivot
	}
	//print_bit_matrix(A, L, L);
	/*
	//printf("remain of inactive rows = %u\n", b);
	if (b != 0){
		printf("Cannot solve equation, un-resolved = %u\n", b);
		free(order_list);
		free(A);
		return 9;
	}
	*/

	// (d) Back-substitution
	for (row_id = 0; row_id < L; row_id++){
		//printf("degree_list[%u] = %u, order_list = %u\n", i, degree_list[i], order_list[i]);
		if (degree_list[row_id] != 0xfffd)	// skip old pivot rows
			continue;

		// XOR solved row to old pivot rows
		col_id = order_list[row_id];	// load order of rows
		C_p = Cd + (size_t)Tal * row_id;
		row_p = A;
		for (i = 0; i < L; i++){
			if (degree_list[i] == 0xffff){
				if (bitmask_check(row_p, col_id) != 0)	// XOR symbols
					align_xor(Cd + (size_t)Tal * i, C_p, Tal);
			}
			row_p += int_col;
		}
	}

	// Because encoder doesn't fail, no need matrix anymore.
	free(A);

	// restore intermediate symbols
	C_p = rc->C;
	for (row_id = 0; row_id < L; row_id++){
		col_id = order_list[row_id];	// load order of rows
		memcpy(C_p + (size_t)Tal * col_id, Cd, Tal);

		Cd += Tal;	// goto next symbol
	}
	free(order_list);

	return 0;
}

