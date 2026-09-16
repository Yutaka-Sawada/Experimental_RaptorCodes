/*
MIT License
Copyright (c) 2026 Yutaka Sawada
*/
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#include "raptor.h"
#include "store.h"


// returns 0 = success, others = fail
int allocate_symbol_store(raptor_coder *rc)
{
	unsigned int i, SH = rc->S + rc->H;
	unsigned int max_tmp = rc->max_tmp;
	unsigned int *tmp_list;

	rc->tmp_buf = malloc((size_t)rc->Tal * max_tmp);
	if (rc->tmp_buf == NULL)
		return 1;

	// zero fill list at first
	tmp_list = calloc(max_tmp, sizeof(unsigned int));
	if (tmp_list == NULL){
		free(rc->tmp_buf);
		rc->tmp_buf = NULL;
		return 2;
	}
	rc->tmp_list = tmp_list;

	// zero fill the first S+H symbols
	memset(rc->tmp_buf, 0, (size_t)rc->Tal * SH);
	for (i = 0; i < SH; i++)
		tmp_list[i] = rc->L + RAPTOR_EMAX + i;
	rc->num_tmp = SH;

	return 0;
}

void reduce_symbol_store(raptor_coder *rc)
{
	void *tmp_p;
	tmp_p = realloc(rc->tmp_buf, rc->Tal);
	if (tmp_p != NULL){
		rc->tmp_buf = tmp_p;
		rc->max_tmp = 1;
	}
	if (rc->tmp_list){
		free(rc->tmp_list);
		rc->tmp_list = NULL;
	}
	rc->num_tmp = 0;
}

// It returns empty position in temporary buffer.
// 0xffffffff (-1) = fail to enlarge buffer size
unsigned int search_empty_pos(raptor_coder *rc)
{
	unsigned int empty_pos;
	unsigned int last_pos = rc->num_tmp;
	unsigned int max_tmp = rc->max_tmp;
	unsigned int *tmp_list = rc->tmp_list;

	if (last_pos < max_tmp){
		// check the last position at first.
		if (tmp_list[last_pos] == 0)
			return last_pos;

		// search empty position in temporary buffer.
		for (empty_pos = 0; empty_pos < max_tmp; empty_pos++){
			if (tmp_list[empty_pos] == 0)
				return empty_pos;
		}
	}

	// buffer is full.
	void *tmp_p;
	unsigned int add_num = rc->S + rc->H;

	// enlarge temporary buffer
	tmp_p = realloc(rc->tmp_buf, (size_t)rc->Tal * (max_tmp + add_num));
	if (tmp_p == NULL)
		return 0xffffffff;
	rc->tmp_buf = tmp_p;

	// enlarge list of temporary buffer
	tmp_list = realloc(rc->tmp_list, sizeof(unsigned int) * (max_tmp + add_num));
	if (tmp_list == NULL)
		return 0xffffffff;
	rc->tmp_list = tmp_list;
	memset(tmp_list + max_tmp, 0, sizeof(unsigned int) * add_num);	// zero fill new area

	//printf("enlarge temporary store from %u to %u\n", max_tmp, max_tmp + add_num);
	rc->max_tmp = max_tmp + add_num;
	return max_tmp;
}

// unique ID of encoding symbols = L + ESI
// unique ID of check symbols    = L + EMAX ~ L + EMAX + S + H - 1
void put_symbol_store(
	raptor_coder *rc,
	void *data,
	unsigned int store_pos,
	unsigned int uid)		// unique ID
{
	memcpy(rc->tmp_buf + (size_t)rc->Tal * store_pos, data, rc->T);
	rc->tmp_list[store_pos] = uid;
	rc->num_tmp += 1;
}

