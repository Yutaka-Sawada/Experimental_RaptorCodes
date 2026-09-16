#ifndef RAPTOR_STORE_H
#define RAPTOR_STORE_H

int allocate_symbol_store(raptor_coder *rc);

void reduce_symbol_store(raptor_coder *rc);

unsigned int search_empty_pos(raptor_coder *rc);

void put_symbol_store(
	raptor_coder *rc,
	void *data,
	unsigned int store_pos,
	unsigned int uid);		// unique ID

#endif
