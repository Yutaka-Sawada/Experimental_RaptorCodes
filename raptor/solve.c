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
#include "solve.h"
#include "util.h"


/*
Peeling Algorithm
Try to recover intermediate symbols by searching nodes with single edge

This will require several hundred overheads to recover all intermediate symbols.

It returns zero, when it could not recover any symbols.
Or else, it returns non-zero for partial recovery.
*/
unsigned int solve_bg_peel(raptor_coder *rc)
{
	unsigned int i, id;
	unsigned int node_start, node_id, node_pos, node_len;
	unsigned int edge_count, lost_id;
	unsigned int slide_size;
	unsigned int Tal = rc->Tal;
	unsigned char *C = rc->C;	// array of intermediate symbols
	uint64_t *pre_mask = rc->pre_mask;
	unsigned int list_len = rc->num_bg;
	unsigned int *node_list = rc->bg;
	unsigned char *tmp_buf = rc->tmp_buf;
	unsigned int *tmp_list = rc->tmp_list;

	// test every nodes
	slide_size = 0;
	node_start = 0;
	while (node_start < list_len){
		node_id = node_list[node_start];
		node_pos = node_list[node_start + 1];
		node_len = node_list[node_start + 2];

		// test neighbor node
		edge_count = 0;
		for (i = 0; i < node_len; i++){
			id = node_list[node_start + 3 + i];
			if (bitmask_check(pre_mask, id) == 0){	// The neighbor symbol is missing.
				edge_count++;
				lost_id = id;
			} else {
				// If the neighbor symbol exists, XOR symbol and remove the edge.
				align_xor(tmp_buf + (size_t)Tal * node_pos, C + (size_t)Tal * id, Tal);
				//printf("XOR neighbor symbol %u -> stored symbol %u\n", id, node_id);
				node_list[node_start + 3 + i] = ERASE_ID;	// eraser mark
			}
		}
		//printf("node = %u, start = %u, edge = %u / %u\n", node_id, node_start, edge_count, node_len);

		// If there is only one edge, the data is intermediate symbol.
		if (edge_count == 1){
			// copy stored symbol to aligned buffer
			memcpy(C + (size_t)Tal * lost_id, tmp_buf + (size_t)Tal * node_pos, Tal);
			//printf("intermediate symbol %u is recovered from node %u\n", lost_id, node_id);
			bitmask_set(pre_mask, lost_id);
			rc->recover += 1;	// count recovered symbols
			edge_count = 0;
		}

		if (edge_count == 0){	// remove this node from bipartite graph
			tmp_list[node_pos] = 0;	// remove stored symbol
			rc->num_tmp -= 1;
			slide_size += 3 + node_len;
		} else {
			if (slide_size > 0){	// slide this node to front
				node_list[node_start - slide_size] = node_id;
				node_list[node_start - slide_size + 1] = node_pos;
			}
			if ((slide_size > 0) || (edge_count < node_len)){
				node_list[node_start - slide_size + 2] = edge_count;
				id = node_start - slide_size + 3;
				for (i = 0; i < node_len; i++){
					lost_id = node_list[node_start + 3 + i];
					if (lost_id != ERASE_ID){	// skip erased edge
						node_list[id] = lost_id;
						id++;
					}
				}
			}
			if (edge_count < node_len)
				slide_size += node_len - edge_count;
		}
		node_start += 3 + node_len;	// goto next node
	}

	if (slide_size > 0){
		//printf("slide_size = %u\n", slide_size);
		rc->num_bg = list_len - slide_size;
		//print_bg(rc);
	}
	return slide_size;
}

/*
Try to solve by Gaussian Elimination
Gaussian Elimination is very slow. O(n power 3)
Before call this, solve_bg_peel() must return 0 to ensure that left nodes don't exist.

It returns zero for successful recovery.
Or else, it returns non-zero for error.
*/
unsigned int solve_bg_ge(raptor_coder *rc)
{
	unsigned int i, id;
	unsigned int num_col, num_row, num_lost;
	unsigned int node_start, node_id, node_pos, node_len;
	unsigned int L = rc->L;
	unsigned int list_len = rc->num_bg;
	unsigned int *node_list = rc->bg;
	unsigned int *map_col, *map_rev;

	num_lost = L - rc->recover;	// number of lost symbols
	if (num_lost > rc->num_tmp)
		return 0;	// exit suddenly, when missing symbols are more than stored symbols.
	//printf("total = %u, exist = %u, lost = %u\n", L, rc->recover, num_lost);
	//print_bg(rc);

	// map of symbol ID to column ID of matrix
	map_col = calloc(L, sizeof(unsigned int));
	if (map_col == NULL)
		return 1;

	// count number of nodes and their edges
	num_row = 0;
	node_start = 0;
	while (node_start < list_len){
		node_id = node_list[node_start];
		node_len = node_list[node_start + 2];

		// test edges
		for (i = 0; i < node_len; i++){
			id = node_list[node_start + 3 + i];
			map_col[id] += 1;
		}
		num_row++;

		node_start += 3 + node_len;	// goto next node
	}
	if (num_row != rc->num_tmp){	// each row must has its stored symbol
		//printf("num_row = %u, num_tmp = %u\n", num_row, rc->num_tmp);
		free(map_col);
		return 0;	// exit without modifying bipartite graph
	}

/*
	printf("count %u symbol IDs in %u rows -> %u symbols lost\n", L, num_row, num_lost);
	for (i = 0; i < L; i++)
		printf(" %u", map_col[i]);
	printf("\n");
*/

	// reverse map of column ID of matrix to symbol ID
	map_rev = malloc(sizeof(unsigned int) * num_lost);
	if (map_rev == NULL){
		free(map_col);
		return 2;
	}

	// count number of columns
	num_col = 0;
	for (i = 0; i < L; i++){
		id = map_col[i];
		if (id == 0){
			//map_col[i] = 0xffffffff;	// mark of ignored symbol for debug print
			continue;
		}
		// missing symbol
		map_col[i] = num_col;
		map_rev[num_col] = i;
		num_col++;
	}
	//printf("num_row = %u, num_col = %u\n", num_row, num_col);
	if (num_col < num_lost){	// available symbols are too few
		free(map_col);
		free(map_rev);
		return 0;	// exit without modifying bipartite graph
	}

/*
	// mapping lost symbol ID to column ID of matrix
	printf("map %u symbol IDs -> %u column IDs\n", L, num_col);
	for (i = 0; i < L; i++)
		printf(" %d", map_col[i]);
	printf("\n");
	printf("reverse map %u column IDs -> %u symbol IDs\n", num_col, L);
	for (i = 0; i < num_col; i++)
		printf(" %u", map_rev[i]);
	printf("\n");
*/

	// allocate three lists at once
	unsigned int *store_list;
	unsigned short *degree_list, *order_list;
	store_list = malloc(sizeof(unsigned int) * num_row * 2);
	if (store_list == NULL){
		free(map_col);
		free(map_rev);
		return 3;
	}
	degree_list = (unsigned short *)(store_list + num_row);
	order_list = degree_list + num_row;

	// allocate memory for matrix
	unsigned int col_id, row_id;
	unsigned int int_col = (num_col + IDXBITS - 1) / IDXBITS;	// number of integers for columns
	uint64_t *matrix, *row_p;
	matrix = calloc(int_col * num_row, sizeof(uint64_t));
	if (matrix == NULL){
		free(map_col);
		free(map_rev);
		free(store_list);
		return 4;
	}
	//printf("bit matrix's %u columns -> %u integers, %u rows, %g KB\n", num_col, int_col, num_row, (double)(int_col * num_row) / 125.0);

	// construct matrix from bipartite graph
	row_id = 0;
	row_p = matrix;
	node_start = 0;
	while (node_start < list_len){
		node_id = node_list[node_start];
		node_pos = node_list[node_start + 1];
		node_len = node_list[node_start + 2];

		// test node
		degree_list[row_id] = node_len;
		store_list[row_id] = node_pos;	// store position in tmp_buf
		row_id++;

		// test edges
		for (i = 0; i < node_len; i++){
			id = node_list[node_start + 3 + i];
			bitmask_set(row_p, map_col[id]);
		}
		row_p += int_col;

		node_start += 3 + node_len;	// goto next node
	}
	free(map_col);	// No need this list anymore
	//print_bit_matrix(matrix, num_col, num_row);

	// solve equation
	unsigned int pivot_id, pivot_degree;
	unsigned int row_degree, next_id;
	unsigned int Tal = rc->Tal;
	unsigned char *C = rc->C;	// array of intermediate symbols
	uint64_t *pre_mask = rc->pre_mask;
	unsigned char *tmp_buf = rc->tmp_buf;
	unsigned int *tmp_list = rc->tmp_list;
	unsigned char *src_p;
	uint64_t *pivot_p;

	// find a row with the smallest degree
	pivot_degree = 0xffff;
	for (i = 0; i < num_row; i++){
		row_degree = degree_list[i];
		//printf("degree_list[%u] = %u\n", i, row_degree);
		if (pivot_degree > row_degree){
			pivot_degree = row_degree;
			pivot_id = i;
		}
	}

	num_lost = 0;	// count number of resolved rows
	while (pivot_degree < 0xffff){
		col_id = bitmask_ntz(matrix + int_col * pivot_id, int_col);
		//printf("pivot_id = %u, pivot_degree = %u, col_id = %u\n", pivot_id, pivot_degree, col_id);
		pivot_p = matrix + int_col * pivot_id;
		order_list[pivot_id] = col_id;	// save order of rows
		degree_list[pivot_id] = 0xffff;	// no need degree of pivot row

		// XOR pivot_row to other rows
		pivot_degree = 0xffff;
		src_p = tmp_buf + (size_t)Tal * store_list[pivot_id];
		row_p = matrix;
		for (row_id = 0; row_id < num_row; row_id++){
			if (row_id == pivot_id){
				row_p += int_col;
				continue;
			}

			if (bitmask_check(row_p, col_id) != 0){	// XOR rows
				align_xor(tmp_buf + (size_t)Tal * store_list[row_id], src_p, Tal);
				if (degree_list[row_id] != 0xffff){	// re-calculate degree
//					for (i = 0; i < int_col; i++)
//						row_p[i] ^= pivot_p[i];
//					row_degree = bitmask_pop(row_p, int_col);
					row_degree = align_xor_pop(row_p, pivot_p, int_col);
					if (row_degree == 0){	// If row_degree is zero, it won't select this.
						row_degree = 0xffff;
						order_list[row_id] = 0xffff;	// ignore this row at recovery
						//printf("degree_list[%u] = 0, store %u\n", row_id, store_list[row_id]);
						tmp_list[ store_list[row_id] ] = 0;	// remove stored symbol
						rc->num_tmp -= 1;
					}
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
		num_lost++;

/*
		printf("next_id = %u, pivot_degree = %u\n", next_id, pivot_degree);
		if (num_lost >= 11)
			print_bit_matrix(matrix, num_col, num_row);
		if (num_lost >= 12)
			break;
*/

		pivot_id = next_id;	// set next pivot
	}
	//print_bit_matrix(matrix, num_col, num_row);

	//printf("number of resolved rows = %u / %u\n", num_lost, num_col);
	if (num_lost < num_col){
		//printf("Cannot solve equation, un-resolved = %u\n", num_col - num_lost);
		// recover missing symbols if possible
		node_len = 0;	// count number of edges
		num_lost = 0;	// count number of nodes
		row_p = matrix;
		for (row_id = 0; row_id < num_row; row_id++){
			row_degree = bitmask_pop(row_p, int_col);
			//printf("degree_list[%u] = %u\n", row_id, row_degree);
			if (row_degree == 1){
				col_id = order_list[row_id];	// load order of rows
				id = map_rev[col_id];
				node_pos = store_list[row_id];
				memcpy(C + (size_t)Tal * id, tmp_buf + (size_t)Tal * node_pos, Tal);
				//printf("intermediate symbol %u (column %u) is recovered from row %u\n", id, col_id, row_id);
				bitmask_set(pre_mask, id);
				rc->recover += 1;	// count recovered symbols
				tmp_list[node_pos] = 0;	// remove stored symbol
				rc->num_tmp -= 1;
				if (rc->recover == L){
					free(map_rev);
					free(store_list);
					free(matrix);
					return 0;
				}
				row_degree = 0;
			}
			degree_list[row_id] = row_degree;
			node_len += row_degree;
			if (row_degree > 0)
				num_lost++;
			row_p += int_col;
		}
		//printf("edges = %u, nodes = %u\n", node_len, num_lost);
		// re-allocate memory
		list_len = node_len + num_lost * 3;
		list_len += (3 + 40) * 2;	// may add 2 with max degree
		//printf("current max = %u, new max = %u\n", rc->max_bg, list_len);
		if (list_len > rc->max_bg){
			free(rc->bg);
			rc->bg = malloc(list_len * sizeof(unsigned int));
			if (rc->bg == NULL){
				free(map_rev);
				free(store_list);
				free(matrix);
				return 5;
			}
			node_list = rc->bg;
			rc->max_bg = list_len;
			rc->num_bg = 0;
		}
		// re-construct bipartite graph from matrix
		uint64_t tmp;
		node_start = 0;
		row_p = matrix;
		for (row_id = 0; row_id < num_row; row_id++){
			node_len = degree_list[row_id];
			if (node_len == 0){
				row_p += int_col;
				continue;
			}
			node_pos = store_list[row_id];
			node_id = tmp_list[node_pos];
			//printf("row_id = %u, node_id = %u, node_pos = %u, node_len = %u\n", row_id, node_id, node_pos, node_len);
			node_list[node_start] = node_id;
			node_list[node_start + 1] = node_pos;
			node_list[node_start + 2] = node_len;
			node_start += 3;
			for (id = 0; id < int_col; id++){
				tmp = row_p[id];
				if (tmp == 0)
					continue;
				for (i = 0; i < IDXBITS; i++){
					if (tmp & 1){
						col_id = id * IDXBITS + i;
						pivot_id = map_rev[col_id];
						//printf("col_id = %u, sumbol id = %u\n", col_id, pivot_id);
						node_list[node_start] = pivot_id;
						node_start++;
					}
					tmp = tmp >> 1;
				}
			}
			row_p += int_col;
		}
		rc->num_bg = node_start;
		//print_bg(rc);
		// release memory
		free(map_rev);
		free(store_list);
		free(matrix);
		return 0;

	} else {	// If it solved succesfully, no need matrix anymore.
		free(matrix);
	}

	// recover missing intermediate symbols
	for (row_id = 0; row_id < num_row; row_id++){
		col_id = order_list[row_id];	// load order of rows
		if (col_id != 0xffff){
			id = map_rev[col_id];
			node_pos = store_list[row_id];
			memcpy(C + (size_t)Tal * id, tmp_buf + (size_t)Tal * node_pos, Tal);
			//printf("intermediate symbol %u (column %u) is recovered from row %u\n", id, col_id, row_id);
			bitmask_set(pre_mask, id);
			rc->recover += 1;	// count recovered symbols
			tmp_list[node_pos] = 0;	// remove stored symbol
			rc->num_tmp -= 1;
			if (rc->recover == L)
				break;
		}
	}

	free(map_rev);
	free(store_list);
	return 0;
}

/*
refer to a papar by Francisco Lazaro and Gerhard Bauch;
Inactivation Decoding of LT and Raptor Codes: Analysis and Code Design

Try to solve by Inactivation Decoding
Inactivation Decoding is faster than Gaussian Elimination. (30% ~ 40% improvement)
Before call this, solve_bg_peel() must return 0 to ensure that left nodes don't exist.

It returns zero for successful recovery.
Or else, it returns non-zero for error.
*/
unsigned int solve_bg_inact(raptor_coder *rc)
{
	unsigned int i, id;
	unsigned int num_col, num_row, num_lost;
	unsigned int node_start, node_id, node_pos, node_len;
	unsigned int L = rc->L;
	unsigned int list_len = rc->num_bg;
	unsigned int *node_list = rc->bg;
	unsigned int *map_col, *map_rev;

	num_lost = L - rc->recover;	// number of lost symbols
	if (num_lost > rc->num_tmp)
		return 0;	// exit suddenly, when missing symbols are more than stored symbols.
	//printf("total = %u, exist = %u, lost = %u\n", L, rc->recover, num_lost);
	//print_bg(rc);

	// map of symbol ID to column ID of matrix
	map_col = calloc(L, sizeof(unsigned int));
	if (map_col == NULL)
		return 1;

	// count number of nodes and their edges
	num_row = 0;
	node_start = 0;
	while (node_start < list_len){
		node_id = node_list[node_start];
		node_len = node_list[node_start + 2];

		// test edges
		for (i = 0; i < node_len; i++){
			id = node_list[node_start + 3 + i];
			map_col[id] += 1;
		}
		num_row++;

		node_start += 3 + node_len;	// goto next node
	}
	if (num_row != rc->num_tmp){	// each row must has its stored symbol
		//printf("num_row = %u, num_tmp = %u\n", num_row, rc->num_tmp);
		free(map_col);
		return 0;	// exit without modifying bipartite graph
	}

/*
	printf("count %u symbol IDs in %u rows -> %u symbols lost\n", L, num_row, num_lost);
	for (i = 0; i < L; i++)
		printf(" %u", map_col[i]);
	printf("\n");
*/

	// reverse map of column ID of matrix to symbol ID
	map_rev = malloc(sizeof(unsigned int) * num_lost);
	if (map_rev == NULL){
		free(map_col);
		return 2;
	}

	// count number of columns
	num_col = 0;
	for (i = 0; i < L; i++){
		id = map_col[i];
		if (id == 0){
			//map_col[i] = 0xffffffff;	// mark of ignored symbol for debug print
			continue;
		}
		// missing symbol
		map_col[i] = num_col;
		map_rev[num_col] = i;
		num_col++;
	}
	//printf("num_row = %u, num_col = %u\n", num_row, num_col);
	if (num_col < num_lost){	// available symbols are too few
		free(map_col);
		free(map_rev);
		return 0;	// exit without modifying bipartite graph
	}

/*
	// mapping lost symbol ID to column ID of matrix
	printf("map %u symbol IDs -> %u column IDs\n", L, num_col);
	for (i = 0; i < L; i++)
		printf(" %d", map_col[i]);
	printf("\n");
	printf("reverse map %u column IDs -> %u symbol IDs\n", num_col, L);
	for (i = 0; i < num_col; i++)
		printf(" %u", map_rev[i]);
	printf("\n");
*/

	// allocate three lists at once
	unsigned int *store_list;
	unsigned short *degree_list, *order_list;
	store_list = malloc(sizeof(unsigned int) * num_row * 2);
	if (store_list == NULL){
		free(map_col);
		free(map_rev);
		return 3;
	}
	degree_list = (unsigned short *)(store_list + num_row);
	order_list = degree_list + num_row;

	// allocate memory for matrix
	unsigned int col_id, row_id;
	unsigned int int_col = (num_col + IDXBITS - 1) / IDXBITS;	// number of integers for columns
	uint64_t *matrix, *row_p;
	matrix = calloc(int_col * num_row, sizeof(uint64_t));
	if (matrix == NULL){
		free(map_col);
		free(map_rev);
		free(store_list);
		return 4;
	}
	//printf("bit matrix's %u columns -> %u integers, %u rows, %g KB\n", num_col, int_col, num_row, (double)(int_col * num_row) / 125.0);

	// construct matrix from bipartite graph
	row_id = 0;
	row_p = matrix;
	node_start = 0;
	while (node_start < list_len){
		node_id = node_list[node_start];
		node_pos = node_list[node_start + 1];
		node_len = node_list[node_start + 2];

		// test node
		degree_list[row_id] = node_len;
		store_list[row_id] = node_pos;	// store position in tmp_buf
		row_id++;

		// test edges
		for (i = 0; i < node_len; i++){
			id = node_list[node_start + 3 + i];
			bitmask_set(row_p, map_col[id]);
		}
		row_p += int_col;

		node_start += 3 + node_len;	// goto next node
	}
	free(map_col);	// No need this list anymore
	//print_bit_matrix(matrix, num_col, num_row);

	// solve equation
	unsigned int pivot_id, pivot_degree, row_degree;
	unsigned int next_pivot, next_col;
	unsigned int Tal = rc->Tal;
	unsigned char *C = rc->C;	// array of intermediate symbols
	uint64_t *pre_mask = rc->pre_mask;
	unsigned char *tmp_buf = rc->tmp_buf;
	unsigned int *tmp_list = rc->tmp_list;
	unsigned char *src_p;
	uint64_t *act_mask, *pivot_p;

	// setup list of intermediate symbols
	act_mask = (uint64_t *)calloc(int_col, sizeof(uint64_t));
	if (act_mask == NULL){
		free(map_rev);
		free(store_list);
		free(matrix);
		return 5;
	}

	// find a row with the smallest degree
	pivot_degree = 0xfffd;
	for (i = 0; i < num_row; i++){
		row_degree = degree_list[i];
		//printf("degree_list[%u] = %u\n", i, row_degree);
		if (pivot_degree > row_degree){
			pivot_degree = row_degree;
			pivot_id = i;
		}
	}
	col_id = bitmask_ntz(matrix + int_col * pivot_id, int_col);

	// (a) Triangulation process
	// (b) Zero matrix procedure
	num_lost = 0;	// count number of resolvable rows
	while (pivot_degree < 0xfffd){
		/*
		if (col_id >= num_col){
			print_bit_matrix(matrix, num_col, num_row);
			printf("Cannot solve equation, pivot_id = %u, col_id = %u\n", pivot_id, col_id);
			free(map_rev);
			free(store_list);
			free(matrix);
			free(act_mask);
			return 6;
		}
		*/
		//printf("pivot_id = %u, pivot_degree = %u, col_id = %u\n", pivot_id, pivot_degree, col_id);
		pivot_p = matrix + int_col * pivot_id;
		order_list[pivot_id] = col_id;	// save order of rows
		degree_list[pivot_id] = 0xffff;	// no need degree of pivot row

		if (pivot_degree >= 2){	// inactivate following symbols
			for (i = col_id / IDXBITS; i < int_col; i++)
				act_mask[i] |= pivot_p[i];
			//print_bit_matrix(act_mask, num_col, 1);
		}

		// XOR pivot_row to other rows
		pivot_degree = 0xfffd;
		src_p = tmp_buf + (size_t)Tal * store_list[pivot_id];
		row_p = matrix;
		for (row_id = 0; row_id < num_row; row_id++){
			if (degree_list[row_id] >= 0xfffd){	// skip old pivot rows, inactive rows, and empty rows
				row_p += int_col;
				continue;
			}

			if (bitmask_check(row_p, col_id) != 0){	// XOR rows
				align_xor(tmp_buf + (size_t)Tal * store_list[row_id], src_p, Tal);
				// re-calculate degree
				row_degree = align_xor_pop(row_p, pivot_p, int_col);
				if (row_degree == 0){	// If row_degree is zero, it won't select this.
					row_degree = 0xfffe;
					order_list[row_id] = 0xffff;	// ignore this row at recovery
					//printf("degree_list[%u] = 0, store %u\n", row_id, store_list[row_id]);
					tmp_list[ store_list[row_id] ] = 0;	// remove stored symbol
					rc->num_tmp -= 1;
				}
				if (pivot_degree > row_degree){
					// test inactive row
					id = align_test_include(act_mask, row_p, int_col);
					if (id == 0xffffffff){	// all symbols are inactive
						row_degree = 0xfffd;
						//printf("inactivate XORed row %u\n", row_id);
					} else {
						pivot_degree = row_degree;
						next_pivot = row_id;
						next_col = id;
					}
				}
				degree_list[row_id] = row_degree;
				//printf("degree_list[%u] = %u\n", row_id, row_degree);

			} else {	// find a row with the smallest degree
				row_degree = degree_list[row_id];
				if (pivot_degree > row_degree){
					// test inactive row
					id = align_test_include(act_mask, row_p, int_col);
					if (id == 0xffffffff){	// all symbols are inactive
						degree_list[row_id] = 0xfffd;
						//printf("inactivate row %u\n", row_id);
					} else {
						pivot_degree = row_degree;
						next_pivot = row_id;
						next_col = id;
					}
				}
			}
			row_p += int_col;
		}
		num_lost++;

/*
		printf("next_pivot = %u, next_col = %u, pivot_degree = %u\n", next_pivot, next_col, pivot_degree);
		if (num_lost >= 11)
			print_bit_matrix(matrix, num_col, num_row);
		if (num_lost >= 12)
			break;
*/

		pivot_id = next_pivot;	// set next pivot
		col_id = next_col;
	}
	free(act_mask);	// No need to check active symbols anymore
	//print_bit_matrix(matrix, num_col, num_row);

	// find a row with the smallest degree
	pivot_degree = 0xfffd;
	row_p = matrix;
	for (i = 0; i < num_row; i++){
		row_degree = degree_list[i];
		if (row_degree == 0xfffd){	// re-calculate degree
			row_degree = bitmask_pop(row_p, int_col);
			degree_list[i] = row_degree;
		}
		//printf("degree_list[%u] = %u, order_list = %u\n", i, row_degree, order_list[i]);
		if (pivot_degree > row_degree){
			pivot_degree = row_degree;
			pivot_id = i;
		}
		row_p += int_col;
	}

	// (c) Gaussian elimination
	while (pivot_degree < 0xfffd){
		num_lost++;
		col_id = bitmask_ntz(matrix + int_col * pivot_id, int_col);
		/*
		if (col_id >= L){
			printf("Cannot solve equation, pivot_id = %u, col_id = %u\n", pivot_id, col_id);
			free(map_rev);
			free(store_list);
			free(matrix);
			return 7;
		}
		*/
		//printf("pivot_id = %u, pivot_degree = %u, col_id = %u\n", pivot_id, pivot_degree, col_id);
		pivot_p = matrix + int_col * pivot_id;
		order_list[pivot_id] = col_id;	// save order of rows
		degree_list[pivot_id] = 0xfffd;	// pivot row of inactive symbols

		// XOR pivot_row to other rows
		pivot_degree = 0xfffd;
		src_p = tmp_buf + (size_t)Tal * store_list[pivot_id];
		row_p = matrix;
		for (row_id = 0; row_id < num_row; row_id++){
			if ((row_id == pivot_id) || (degree_list[row_id] > 0xfffd)){	// skip old pivot rows and empty rows
				row_p += int_col;
				continue;
			}

			if (bitmask_check(row_p, col_id) != 0){	// XOR rows
				align_xor(tmp_buf + (size_t)Tal * store_list[row_id], src_p, Tal);
				if (degree_list[row_id] != 0xfffd){	// re-calculate degree
					row_degree = align_xor_pop(row_p, pivot_p, int_col);
					if (row_degree == 0){	// If row_degree is zero, it won't select this.
						row_degree = 0xfffe;
						order_list[row_id] = 0xffff;	// ignore this row at recovery
						//printf("degree_list[%u] = 0, store %u\n", row_id, store_list[row_id]);
						tmp_list[ store_list[row_id] ] = 0;	// remove stored symbol
						rc->num_tmp -= 1;
					}
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

		pivot_id = next_pivot;	// set next pivot
	}
	//print_bit_matrix(matrix, num_col, num_row);

	//printf("number of resolvable rows = %u / %u\n", num_lost, num_col);
	if (num_lost < num_col){
		//printf("Cannot solve equation, un-resolved = %u\n", num_col - num_lost);
		// recover missing symbols if possible
		node_len = 0;	// count number of edges
		num_lost = 0;	// count number of nodes
		row_p = matrix;
		for (row_id = 0; row_id < num_row; row_id++){
			row_degree = bitmask_pop(row_p, int_col);
			//printf("degree_list[%u] = %u\n", row_id, row_degree);
			if (row_degree == 1){
				col_id = order_list[row_id];	// load order of rows
				id = map_rev[col_id];
				node_pos = store_list[row_id];
				memcpy(C + (size_t)Tal * id, tmp_buf + (size_t)Tal * node_pos, Tal);
				//printf("intermediate symbol %u (column %u) is recovered from row %u\n", id, col_id, row_id);
				bitmask_set(pre_mask, id);
				rc->recover += 1;	// count recovered symbols
				tmp_list[node_pos] = 0;	// remove stored symbol
				rc->num_tmp -= 1;
				if (rc->recover == L){
					free(map_rev);
					free(store_list);
					free(matrix);
					return 0;
				}
				row_degree = 0;
			}
			degree_list[row_id] = row_degree;
			node_len += row_degree;
			if (row_degree > 0)
				num_lost++;
			row_p += int_col;
		}
		//printf("edges = %u, nodes = %u\n", node_len, num_lost);
		// re-allocate memory
		list_len = node_len + num_lost * 3;
		list_len += (3 + 40) * 2;	// may add 2 with max degree
		//printf("current max = %u, new max = %u\n", rc->max_bg, list_len);
		if (list_len > rc->max_bg){
			free(rc->bg);
			rc->bg = malloc(list_len * sizeof(unsigned int));
			if (rc->bg == NULL){
				free(map_rev);
				free(store_list);
				free(matrix);
				return 8;
			}
			node_list = rc->bg;
			rc->max_bg = list_len;
			rc->num_bg = 0;
		}
		// re-construct bipartite graph from matrix
		uint64_t tmp;
		node_start = 0;
		row_p = matrix;
		for (row_id = 0; row_id < num_row; row_id++){
			node_len = degree_list[row_id];
			if (node_len == 0){
				row_p += int_col;
				continue;
			}
			node_pos = store_list[row_id];
			node_id = tmp_list[node_pos];
			//printf("row_id = %u, node_id = %u, node_pos = %u, node_len = %u\n", row_id, node_id, node_pos, node_len);
			node_list[node_start] = node_id;
			node_list[node_start + 1] = node_pos;
			node_list[node_start + 2] = node_len;
			node_start += 3;
			for (id = 0; id < int_col; id++){
				tmp = row_p[id];
				if (tmp == 0)
					continue;
				for (i = 0; i < IDXBITS; i++){
					if (tmp & 1){
						col_id = id * IDXBITS + i;
						pivot_id = map_rev[col_id];
						//printf("col_id = %u, sumbol id = %u\n", col_id, pivot_id);
						node_list[node_start] = pivot_id;
						node_start++;
					}
					tmp = tmp >> 1;
				}
			}
			row_p += int_col;
		}
		rc->num_bg = node_start;
		//print_bg(rc);
		// release memory
		free(map_rev);
		free(store_list);
		free(matrix);
		return 0;
	}

	// (d) Back-substitution
	for (row_id = 0; row_id < num_row; row_id++){
		//printf("degree_list[%u] = %u, order_list = %u\n", row_id, degree_list[row_id], order_list[row_id]);
		if (degree_list[row_id] != 0xfffd)	// skip old pivot rows
			continue;

		// XOR solved row to old pivot rows
		col_id = order_list[row_id];	// load order of rows
		src_p = tmp_buf + (size_t)Tal * store_list[row_id];
		row_p = matrix;
		for (i = 0; i < num_row; i++){
			if (degree_list[i] == 0xffff){
				if (bitmask_check(row_p, col_id) != 0)	// XOR symbols
					align_xor(tmp_buf + (size_t)Tal * store_list[i], src_p, Tal);
			}
			row_p += int_col;
		}
	}

	// If it solved succesfully, no need matrix anymore.
	free(matrix);

	// recover missing intermediate symbols
	for (row_id = 0; row_id < num_row; row_id++){
		col_id = order_list[row_id];	// load order of rows
		if (col_id != 0xffff){
			id = map_rev[col_id];
			node_pos = store_list[row_id];
			memcpy(C + (size_t)Tal * id, tmp_buf + (size_t)Tal * node_pos, Tal);
			//printf("intermediate symbol %u (column %u) is recovered from row %u\n", id, col_id, row_id);
			bitmask_set(pre_mask, id);
			rc->recover += 1;	// count recovered symbols
			tmp_list[node_pos] = 0;	// remove stored symbol
			rc->num_tmp -= 1;
			if (rc->recover == L)
				break;
		}
	}

	free(map_rev);
	free(store_list);
	return 0;
}

