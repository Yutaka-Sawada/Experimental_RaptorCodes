#ifndef RAPTOR_BIPARTITE_H
#define RAPTOR_BIPARTITE_H

// mark of erased node or edge
#define ERASE_ID 0xffffffff


int make_bipartite_decode(raptor_coder *rc);

void free_bipartite(raptor_coder *rc);

int bg_add_node(
	raptor_coder *rc,
	unsigned int esi,			// encoding symbol ID
	unsigned int store_pos);	// position in temporary buffer

int bg_add_decode(raptor_coder *rc);

void print_bg(raptor_coder *rc);

#endif
