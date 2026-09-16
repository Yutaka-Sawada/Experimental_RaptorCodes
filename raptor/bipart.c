/*
MIT License
Copyright (c) 2026 Yutaka Sawada
*/
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#include "raptor.h"
#include "bipart.h"
#include "bitmask.h"
#include "gene.h"
#include "precode.h"
#include "store.h"
#include "util.h"


/*
L = K + S + H
unique ID of left node  : intermediate symbols = 0 ~ L - 1
unique ID of right node : encoding symbols = L + ESI
                          check symbols    = L + EMAX ~ L + EMAX + S + H - 1

construction of Bipartite Graph:
{
	unique ID of right node,
	position in temporary buffer,
	length of following items (number of valid IDs may be fewer than this value.),
	unique ID of left node 1,
	unique ID of left node 2,
	...
	unique ID of left node N
} for each node continuously
*/

// returns 0 = success, others = fail
int make_bipartite_encode(raptor_coder *rc)
{
	unsigned int K = rc->K;
	unsigned int S = rc->S;
	unsigned int H = rc->H;
	unsigned int L = rc->L;
	unsigned char *A, *row_p;

	// allocate matrix of L * (S + H)
	A = calloc(L * (S + H), 1);
	if (A == NULL)
		return 1;

	// make constraint matrix of the precode
	unsigned int count;	// number of 1s in matrix
	count = init_precode_matrix(rc, A);
	if (count == 0){
		free(A);
		return 2;
	}
	//printf("number of 1s in constraint matrix = %u\n", count);
	//print_matrix(A, L, S + H);

	// number of check symbols = S + H
	unsigned int num_item = (S + H) * 3 + count;

	// possible edges of source symbols
	num_item += K * (3 + 5) + 40;	// average degree would be 5. spike is 40.
	num_item = (num_item + 15) & ~15;	// multiple of 64-bytes

	// allocate memory for bipartite graph
	rc->bg = malloc(num_item * sizeof(unsigned int));
	if (rc->bg == NULL){
		free(A);
		return 3;
	}
	rc->max_bg = num_item;
	rc->num_bg = 0;

	// copy check symbols from matrix A to bipartite graph
	unsigned int col_id, row_id;
	unsigned int *node_list = rc->bg;
	unsigned int node_start, node_id, node_pos, node_len;

	row_p = A;
	node_start = 0;
	for (row_id = 0; row_id < S + H; row_id++){
		node_id = L + RAPTOR_EMAX + row_id;	// unique ID of check symbol
		node_pos = row_id;					// the first S+H symbols in tmp_buf are check symbols
		node_len = 0;
		for (col_id = 0; col_id < L; col_id++){
			if (row_p[col_id] != 0){
				node_list[node_start + 3 + node_len] = col_id;
				node_len++;
			}
		}
		node_list[node_start    ] = node_id;
		node_list[node_start + 1] = node_pos;
		node_list[node_start + 2] = node_len;
		node_start += 3 + node_len;
		row_p += L;
	}
	rc->num_bg = node_start;
	//print_bg(rc);

	free(A);
	return 0;
}

// returns 0 = success, others = fail
int make_bipartite_decode(raptor_coder *rc)
{
	unsigned int K = rc->K;
	unsigned int S = rc->S;
	unsigned int H = rc->H;
	unsigned int L = rc->L;
	unsigned int N = rc->loaded;	// the number of received encoding symbols
	unsigned char *A, *row_p;

	// allocate matrix of L * (S + H)
	A = calloc(L * (S + H), 1);
	if (A == NULL)
		return 1;

	// make constraint matrix of the precode
	unsigned int count;	// number of 1s in matrix
	count = init_precode_matrix(rc, A);
	if (count == 0){
		free(A);
		return 2;
	}
	//printf("number of 1s in constraint matrix = %u\n", count);
	//print_matrix(A, L, S + H);

	// number of check symbols = S + H
	unsigned int num_item = (S + H) * 3 + count;

	// possible edges of encoding symbols
	num_item += N * (3 + 5) + 40;	// average degree would be 5. spike is 40.
	num_item = (num_item + 15) & ~15;	// multiple of 64-bytes

	// allocate memory for bipartite graph
	rc->bg = malloc(num_item * sizeof(unsigned int));
	if (rc->bg == NULL){
		free(A);
		return 3;
	}
	rc->max_bg = num_item;
	rc->num_bg = 0;

	// copy check symbols from matrix A to bipartite graph
	unsigned int col_id, row_id;
	unsigned int *node_list = rc->bg;
	unsigned int node_start, node_id, node_pos, node_len;

	row_p = A;
	node_start = 0;
	for (row_id = 0; row_id < S + H; row_id++){
		node_id = L + RAPTOR_EMAX + row_id;	// unique ID of check symbol
		node_pos = row_id;					// the first S+H symbols in tmp_buf are check symbols
		node_len = 0;
		for (col_id = 0; col_id < L; col_id++){
			if (row_p[col_id] != 0){
				node_list[node_start + 3 + node_len] = col_id;
				node_len++;
			}
		}
		node_list[node_start    ] = node_id;
		node_list[node_start + 1] = node_pos;
		node_list[node_start + 2] = node_len;
		node_start += 3 + node_len;
		row_p += L;
	}
	rc->num_bg = node_start;
	//print_bg(rc);

	free(A);
	return 0;
}

void free_bipartite(raptor_coder *rc)
{
	if (rc->bg){
		free(rc->bg);
		rc->bg = NULL;
	}
	rc->num_bg = 0;
	rc->max_bg = 0;
}

// add edges of the encoding symbol in bipartite graph
// returns 0 = success, others = fail
int bg_add_node(
	raptor_coder *rc,
	unsigned int esi,			// encoding symbol ID
	unsigned int store_pos)		// position in temporary buffer
{
	unsigned int j;
	unsigned int node_start = rc->num_bg;

	// LT encoding
	unsigned int d, a, b;
	unsigned int K = rc->K;
	unsigned int L = rc->L;
	unsigned int Ld = rc->Ld;

	Trip(K, Ld, esi, &d, &a, &b);
	//printf("esi = %u, d = %u, a = %u, b = %u\n", esi, d, a, b);

	if (node_start + 3 + d >= rc->max_bg){	// buffer is full.
		void *tmp_p;
		// this one and next node with max degree
		unsigned int add_num = 3 + d + 3 + 40;
		add_num = (add_num + 15) & ~15;	// multiple of 64-bytes
		tmp_p = realloc(rc->bg, (size_t)(rc->max_bg + add_num) * sizeof(unsigned int));
		if (tmp_p == NULL)
			return -1;
		//printf("enlarge Bipartite Graph from %u to %u\n", rc->max_bg, rc->max_bg + add_num);
		rc->bg = tmp_p;
		rc->max_bg += add_num;
	}

	// add node
	unsigned int *node_list = rc->bg + node_start;
	unsigned int node_len = 0;
	node_list[0] = L + esi;
	node_list[1] = store_pos;
	node_list[2] = d;	// degree = number of edges

	if (d > L)
		d = L;
	while (b >= L)
		b = (b + a) % Ld;
	node_list[3] = b;
	node_len++;
	for (j = 1; j < d; j++){
		b = (b + a) % Ld;
		while (b >= L)
			b = (b + a) % Ld;
		node_list[3 + node_len] = b;
		node_len++;
	}

	// sort edges
	//qsort(node_list + 3, node_len, sizeof(unsigned int), compare_uint32);
	selection_sort(node_list + 3, node_len);

	rc->num_bg += 3 + node_len;
	return 0;
}

// add edges of received encoding symbols into bipartite graph
// call this once after make_bipartite_decode()
int bg_add_decode(raptor_coder *rc)
{
	unsigned int i;
	unsigned int L = rc->L;
	unsigned int N = rc->loaded;	// the number of received encoding symbols
	unsigned int *esi_list = rc->tmp_list;	// list of ESI
	unsigned int esi, store_pos;

	store_pos = rc->S + rc->H;	// skip check symbols
	for (i = 0; i < N; i++){
		esi = esi_list[store_pos];
/*
No need additional checking ?
		while ((esi == 0) || (esi >= L + RAPTOR_EMAX)){	// ignore check symbols
			store_pos++;
			esi = esi_list[store_pos];
		}
*/
		esi -= L;	// convert from unique ID to ESI
		if (bg_add_node(rc, esi, store_pos) != 0)
			return -1;
		store_pos++;
	}

	//print_bg(rc);
	return 0;
}

void print_bg(raptor_coder *rc)
{
	unsigned int i, j, id;
	unsigned int node_id, node_pos, node_len, edge_count;
	unsigned int max = rc->num_bg;
	unsigned int *list = rc->bg;
	uint64_t *pre_mask = rc->pre_mask;

	if (max == 0)
		return;

	printf(" %u items / max %u\n", max, rc->max_bg);
	i = 0;
	while (i < max){
		node_id = list[i];
		node_pos = list[i + 1];
		node_len = list[i + 2];
//		if (node_id == ERASE_ID){
//			i += 3 + node_len;
//			continue;	// ignore erased node
//		}
		printf("{");
		edge_count = 0;
		for (j = 0; j < node_len; j++){
			id = list[i + 3 + j];
//			if (id == ERASE_ID)
//				continue;	// ignore erased edge
			if (edge_count != 0)
				printf(" ");
			if (bitmask_check(pre_mask, id) == 0){
				printf("%u", id);	// missing
			} else {
				printf("e%u", id);	// exist
			}
			edge_count++;
		}
		printf("} = %u, pos %u\n", node_id, node_pos);
		i += 3 + node_len;
	}
}
