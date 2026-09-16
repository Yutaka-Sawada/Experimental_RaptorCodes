// Raptor Codes
// Copyright (c) 2026 Yutaka Sawada
// MIT License

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#include "raptor.h"
#include "precode.h"
#include "bipart.h"
#include "bitmask.h"
#include "gene.h"
#include "store.h"
#include "solve.h"
#include "util.h"


raptor_coder *raptor_encoder_new(
	unsigned int K,		// number of source symbols
	unsigned int T)		// size of each symbol in bytes
{
	if ((K > RAPTOR_KMAX) || (K < RAPTOR_KMIN) || (T <= 1))
		return NULL;

	raptor_coder *rc = calloc(1, sizeof(raptor_coder));
	if (rc == NULL)
		return NULL;

	rc->K = K;
	rc->T = T;
	rc->Tal = (T + SYMBOL_ALIGN - 1) & ~(SYMBOL_ALIGN - 1);
	init_precode(rc);

	// allocate L intermediate symbols
	size_t Tal = rc->Tal;
	rc->C = malloc(Tal * rc->L);
	if (rc->C == NULL){
		raptor_free(rc);
		return NULL;
	}

	// allocate buffer for K source symbols (and S+H temporary symbols)
	rc->tmp_buf = malloc(Tal * rc->L);
	if (rc->tmp_buf == NULL){
		raptor_free(rc);
		return NULL;
	}
	//rc->max_tmp = rc->L;
	memset(rc->tmp_buf, 0, Tal * (rc->S + rc->H));	// zero fill the first S+H symbols

	// setup list of source symbols
	rc->id_mask = calloc(bitmask_len(K), sizeof(uint64_t));	// zero fill
	if (rc->id_mask == NULL){
		raptor_free(rc);
		return NULL;
	}

	return rc;
}


void raptor_free(raptor_coder *rc)
{
	if (rc == NULL)
		return;
	if (rc->C)
		free(rc->C);
	if (rc->id_mask)
		free(rc->id_mask);
	if (rc->tmp_buf)
		free(rc->tmp_buf);
	if (rc->tmp_list)
		free(rc->tmp_list);
	if (rc->bg)
		free(rc->bg);
	free(rc);
	rc = NULL;
}

// add some source symbols
int raptor_encoder_add(
	raptor_coder *rc,
	void *data,				// input bytes of some adding source symbols
	unsigned int off_id,	// offset to the first Encoding Symbol ID
	unsigned int num_sym)	// number of adding source symbols
{
	unsigned int K = rc->K;
	if (rc->loaded >= K)
		return RAPTOR_ADD_END;	// no need to add more.
	if ((off_id >= K) || (off_id + num_sym > K))
		return RAPTOR_ADD_PARA;	// given block ID or range is invalid.
	if (rc->max_esi != 0)
		return RAPTOR_ADD_DIF;	// decoder is different.

	int ret;
	unsigned int esi;	// Encoding Symbol ID
	unsigned int T = rc->T;
	unsigned int Tal = rc->Tal;
	uint64_t *id_mask = rc->id_mask;
	unsigned char *data_p = data;
	unsigned char *buf_p = rc->tmp_buf + (size_t)Tal * (rc->S + rc->H + off_id);

	ret = RAPTOR_ADD_OK;
	for (esi = off_id; esi < off_id + num_sym; esi++){
		if (bitmask_check(id_mask, esi) != 0){
			ret = RAPTOR_ADD_IGN;	// This source symbol was added already.

		} else {
			// copy source block to aligned buffer
			memcpy(buf_p, data_p, T);
			bitmask_set(id_mask, esi);
			rc->loaded += 1;	// count added symbols
		}

		// goto next source symbol
		buf_p += Tal;
		data_p += T;
	}

	// When all source symbols were added, it creates intermediate symbols.
	if (rc->loaded >= K){
		//printf("added %u symbols OK\n", rc->loaded);

		// make generator matrix of precode and calculate intermediate symbols
		// While Inactivation Decoding is 30~40% faster, it has patent issue.
		if (make_generator_matrix(rc) != 0)		// Gaussian Elimination
		//if (make_inactive_matrix(rc) != 0)	// Inactivation Decoding
			return RAPTOR_ENC_ERR;

		// Reduce temporary buffer for work space
		reduce_symbol_store(rc);
	}

	return ret;
}

// create some encoding symbols at once
// When esi < number of source symbols, source symbol is restored.
// When esi >= number of source symbols, encoding symbol is created.
int raptor_encoder_create(
	raptor_coder *rc,
	void *data,				// output bytes of some created symbols
	unsigned int off_id,	// offset to the first encoding symbol ID
	unsigned int num_sym)	// number of creating symbols
{
	unsigned int K = rc->K;
	if (rc->loaded < K)
		return RAPTOR_ENC_LACK;	// need to add source symbols at first
	if (off_id + num_sym > RAPTOR_EMAX)
		return RAPTOR_ENC_PARA;	// given symbol ID or range is invalid.
	if (rc->max_esi != 0)
		return RAPTOR_ENC_DIF;	// decoder is different.

	unsigned int esi;	// encoding symbol ID
	unsigned int d, a, b, j;
	unsigned int T = rc->T;
	unsigned int Tal = rc->Tal;
	unsigned int L = rc->L;
	unsigned int Ld = rc->Ld;
	unsigned char *data_p = data;
	unsigned char *C = rc->C;
	unsigned char *tmp_buf = rc->tmp_buf;

	// LT encoding
	//unsigned int ave = 0;
	for (esi = off_id; esi < off_id + num_sym; esi++){
		Trip(K, Ld, esi, &d, &a, &b);
		//printf("ESI = %u: d = %u, a = %u, b = %u\n", esi, d, a, b);
		//ave += d;

		if (d > L)
			d = L;
		while (b >= L)
			b = (b + a) % Ld;
		memcpy(tmp_buf, C + (size_t)Tal * b, Tal);	// result = C[b]
		for (j = 1; j < d; j++){
			b = (b + a) % Ld;
			while (b >= L)
				b = (b + a) % Ld;
			align_xor(tmp_buf, C + (size_t)Tal * b, Tal);	// result = result ^ C[b]
		}

		memcpy(data_p, tmp_buf, T);
		data_p += T;	// goto next symbol
	}
	//printf("average degree = %g\n", (double)ave / num_sym);

	return RAPTOR_ENC_OK;
}


raptor_coder *raptor_decoder_new(
	unsigned int K,		// number of source symbols
	unsigned int T,		// size of each symbol in bytes
	unsigned int max_esi)
{
	if ((K > RAPTOR_KMAX) || (K < RAPTOR_KMIN) || (T <= 1))
		return NULL;

	raptor_coder *rc = calloc(1, sizeof(raptor_coder));
	if (rc == NULL)
		return NULL;

	rc->K = K;
	rc->T = T;
	rc->Tal = (T + SYMBOL_ALIGN - 1) & ~(SYMBOL_ALIGN - 1);
	init_precode(rc);

	// allocate L intermediate symbols
	rc->C = malloc((size_t)rc->Tal * rc->L);
	if (rc->C == NULL){
		raptor_free(rc);
		return NULL;
	}

	// If max_esi is strange, set default value automatically.
	if ((max_esi >= RAPTOR_EMAX || max_esi < K)){
		rc->max_esi = K * 2 + RAPTOR_OVERHEAD - 1;
		rc->max_tmp = rc->L + RAPTOR_OVERHEAD;
	} else {
		rc->max_esi = max_esi;
		rc->max_tmp = rc->S + rc->H + max_esi + 1;
		if (rc->max_tmp > rc->L + K)
			rc->max_tmp = rc->L + K;	// 100% overhead should be enough.
	}

	// allocate buffer for temporary data
	if (allocate_symbol_store(rc) != 0){
		raptor_free(rc);
		return NULL;
	}

	// setup list of encoding symbols
	rc->id_mask = calloc(bitmask_len(rc->max_esi + 1), sizeof(uint64_t));
	if (rc->id_mask == NULL){
		raptor_free(rc);
		return NULL;
	}

	// setup list of pre-coding symbols
	rc->pre_mask = calloc(bitmask_len(rc->L), sizeof(uint64_t));
	if (rc->pre_mask == NULL){
		raptor_free(rc);
		return NULL;
	}

	return rc;
}

// Return the number of possible extra symbols to recover lost symbols.
// The number of source symbols and additional symbols will be required.
// The numbers of created repair symbols and lost symbols won't affect.
unsigned int raptor_decoder_overhead(raptor_coder *rc)
{
	// It may require 1.6% overhead at least 2, and at most RAPTOR_OVERHEAD.
	unsigned int num_overhead = (rc->K + 63) / 64;
	if (num_overhead < 2)
		num_overhead = 2;
	if (num_overhead > RAPTOR_OVERHEAD)
		num_overhead = RAPTOR_OVERHEAD;

	return num_overhead;
}

int raptor_decoder_add(
	raptor_coder *rc,
	void *data,				// input bytes of an adding symbol
	unsigned int esi)		// encoding symbol ID
{
	if (rc->max_esi == 0)
		return RAPTOR_ADD_DIF;	// encoder is different.
	if (esi > rc->max_esi)
		return RAPTOR_ADD_PARA;	// given symbol ID is invalid.
	if (rc->recover >= rc->L)
		return RAPTOR_ADD_END;	// no need to add more.

	uint64_t *id_mask = rc->id_mask;
	if (bitmask_check(id_mask, esi) != 0)
		return RAPTOR_ADD_IGN;	// This encoding symbol was added already.

	// search empty position in temporary buffer.
	unsigned int store_pos = search_empty_pos(rc);
	//printf("store_pos = %u / %u\n", store_pos, rc->max_tmp);
	if (store_pos >= rc->max_tmp)
		return RAPTOR_ADD_ERR;	// failed to enlarge temporary buffer.

	// store this encoding symbol in temporary buffer
	memcpy(rc->tmp_buf + (size_t)rc->Tal * store_pos, data, rc->T);
	rc->tmp_list[store_pos] = rc->L + esi;	// unique ID of encoding symbol
	rc->num_tmp += 1;

	if (rc->bg != NULL){	// When graph was made already,
		// add new node into bipartite graph
		if (bg_add_node(rc, esi, store_pos) != 0)
			return RAPTOR_ADD_ERR;
		//print_bg(rc);
	}

	bitmask_set(id_mask, esi);
	rc->loaded += 1;	// count added symbols

	return RAPTOR_ADD_OK;
}

// Try to solve by Iterative Decoding
int raptor_decoder_solve_iterative(raptor_coder *rc)
{
	unsigned int L = rc->L;
	if (rc->recover >= L)
		return RAPTOR_ENC_OK;	// solved already
	if (rc->loaded < rc->K)
		return RAPTOR_ENC_LACK;	// need more symbols
	if (rc->max_esi == 0)
		return RAPTOR_ENC_DIF;	// encoder is different.

	if (rc->bg == NULL){	// make bipartite graph at first
		if (make_bipartite_decode(rc) != 0)
			return RAPTOR_ENC_ERR;
		if (bg_add_decode(rc) != 0)
			return RAPTOR_ENC_ERR;
	}
	//print_bg(rc);

	// Iterative Decoding with Peeling Algorithm
	unsigned int prev = rc->recover;
	do {
		if (solve_bg_peel(rc) == 0)
			break;
	} while (rc->recover < L);
	printf("Peeling Decoder recovered %u after %u input, total %u/%u.\n", rc->recover - prev, rc->loaded, rc->recover, L);

	if (rc->recover < L)
		return RAPTOR_ENC_LACK;	// need to add more symbols

	// Erase graph
	free_bipartite(rc);

	// No need to add encoding symbols.
	free(rc->id_mask);
	rc->id_mask = NULL;

	// No need to check intermediate symbols.
	free(rc->pre_mask);
	rc->pre_mask = NULL;

	// Reduce temporary buffer for work space
	reduce_symbol_store(rc);

	return RAPTOR_ENC_OK;
}

// Try to solve by Hybrid Decoding
// This starts faster decoder at first, then finishes with slower decoder next.
int raptor_decoder_solve_hybrid(raptor_coder *rc)
{
	unsigned int L = rc->L;
	if (rc->recover >= L)
		return RAPTOR_ENC_OK;	// solved already
	if (rc->loaded < rc->K)
		return RAPTOR_ENC_LACK;	// need more symbols
	if (rc->max_esi == 0)
		return RAPTOR_ENC_DIF;	// encoder is different.

	if (rc->bg == NULL){	// make bipartite graph at first
		if (make_bipartite_decode(rc) != 0)
			return RAPTOR_ENC_ERR;
		if (bg_add_decode(rc) != 0)
			return RAPTOR_ENC_ERR;
	}
	//print_bg(rc);

	// Iterative Decoding with Peeling Algorithm at first
	unsigned int prev = rc->recover;
	do {
		if (solve_bg_peel(rc) == 0)
			break;
	} while (rc->recover < L);
	printf("Peeling Decoder recovered %u after %u input, total %u/%u.\n", rc->recover - prev, rc->loaded, rc->recover, L);

	// Gaussian Elimination and Inactivation Decoding are alternative each other.
	// Gaussian Elimination is slower.
	// Inactivation Decoding has patent issue.
	if (rc->recover < L){	// Gaussian Elimination once at last
		prev = rc->recover;
		if (solve_bg_ge(rc) != 0)
			return RAPTOR_ENC_ERR;	// lack of memory
		printf("Gaussian Elimination recovered %u after %u input, total %u/%u\n", rc->recover - prev, rc->loaded, rc->recover, L);
	}
/*
	if (rc->recover < L){	// Inactivation Decoding once at last
		prev = rc->recover;
		if (solve_bg_inact(rc) != 0)
			return RAPTOR_ENC_ERR;	// lack of memory
		printf("Inactivation Decoding recovered %d after %u input, total %u/%u.\n", rc->recover - prev, rc->loaded, rc->recover, L);
	}
*/

	if (rc->recover < L)
		return RAPTOR_ENC_LACK;	// need to add more symbols

	// Erase graph
	free_bipartite(rc);

	// No need to add encoding symbols.
	free(rc->id_mask);
	rc->id_mask = NULL;

	// No need to check intermediate symbols.
	free(rc->pre_mask);
	rc->pre_mask = NULL;

	// Reduce temporary buffer for work space
	reduce_symbol_store(rc);

	return RAPTOR_ENC_OK;
}

// recover some source symbols at once
int raptor_decoder_recover(
	raptor_coder *rc,
	void *data,				// output bytes of some recovered symbols
	unsigned int off_id,	// offset to the first encoding symbol ID
	unsigned int num_sym)	// number of recovering symbols
{
	if (rc->max_esi == 0)
		return RAPTOR_ENC_DIF;	// encoder is different.
	unsigned int L = rc->L;
	if (rc->recover < rc->L)
		return RAPTOR_ENC_LACK;	// need to restore intermediate symbols at first
	if (off_id + num_sym > RAPTOR_EMAX)
		return RAPTOR_ENC_PARA;	// given symbol ID or range is invalid.

	unsigned int esi;	// encoding symbol ID
	unsigned int d, a, b, j;
	unsigned int T = rc->T;
	unsigned int Tal = rc->Tal;
	unsigned int K = rc->K;
	unsigned int Ld = rc->Ld;
	unsigned char *data_p = data;
	unsigned char *C = rc->C;
	unsigned char *tmp_buf = rc->tmp_buf;

	// LT encoding
	for (esi = off_id; esi < off_id + num_sym; esi++){
		Trip(K, Ld, esi, &d, &a, &b);
		//printf("ESI = %u: d = %u, a = %u, b = %u\n", esi, d, a, b);

		if (d > L)
			d = L;
		while (b >= L)
			b = (b + a) % Ld;
		memcpy(tmp_buf, C + (size_t)Tal * b, Tal);	// result = C[b]
		for (j = 1; j < d; j++){
			b = (b + a) % Ld;
			while (b >= L)
				b = (b + a) % Ld;
			align_xor(tmp_buf, C + (size_t)Tal * b, Tal);	// result = result ^ C[b]
		}

		memcpy(data_p, tmp_buf, T);
		data_p += T;	// goto next symbol
	}

	return RAPTOR_ENC_OK;
}

