#ifndef RAPTOR_CODES_H
#define RAPTOR_CODES_H

#include <stdint.h>


#define SYMBOL_ALIGN	8	// This value must be multiple of 8.

#define RAPTOR_KMAX	8192	// the maximum number of source symbols
#define RAPTOR_KMIN	4		// a minimum target on the number of symbols
#define RAPTOR_EMAX	65536	// the maximum number of encoding symbols

#define RAPTOR_OVERHEAD	8	// 8 overhead may be enough.

#define RAPTOR_ADD_END   2	// It got enough symbols already. Repair is possible.
#define RAPTOR_ADD_IGN   1	// It ignored the symbol, as it was added already or useless.
#define RAPTOR_ADD_OK    0
#define RAPTOR_ADD_ERR  -1	// Fatal error (mostly insufficient memory)
#define RAPTOR_ADD_DIF  -2	// Either calling encoder or decoder is exclusive.
#define RAPTOR_ADD_PARA -3	// Given parameter is bad.

#define RAPTOR_ENC_LACK  3	// It needs to add more symbols.
#define RAPTOR_ENC_OK    0
#define RAPTOR_ENC_ERR  -1	// Fatal error (mostly insufficient memory)
#define RAPTOR_ENC_DIF  -2	// Either calling encoder or decoder is exclusive.
#define RAPTOR_ENC_PARA -3	// Given parameter is bad.


typedef struct {
	// source symbols
	unsigned int K;			// number of source symbols
	unsigned int T;			// size of each input or output symbol in bytes
	uint64_t *id_mask;		// list of encoding symbol state

	// intermediate symbols
	unsigned int Tal;		// aligned symbol size, which MUST be a multiple of SYMBOL_ALIGN
	unsigned int L;			// number of pre-coding symbols (L = K + S + H)
	unsigned int S;			// number of LDPC symbols
	unsigned int H;			// number of Half symbols
	unsigned int Ld;		// the smallest prime integer greater than or equal to L
	unsigned char *C;		// array of intermediate symbols (0 ~ L-1)
	uint64_t *pre_mask;		// list of pre-coding symbol state

	// work space to store symbols temporary
	unsigned int max_esi;	// maximum ESI of possible encoding symbols
	unsigned char *tmp_buf;	// address of temporary buffer
	unsigned int num_tmp;	// number of stored symbols
	unsigned int max_tmp;	// maximum number of symbols in temporary buffer
	unsigned int *tmp_list;	// list of symbols in temporary buffer

	// nodes of neighbor symbols in bipartite graph
	unsigned int *bg;		// bipartite graph
	unsigned int num_bg;	// number of items
	unsigned int max_bg;	// maximum number of items

	// process counter
	unsigned int loaded;	// number of added symbols
	unsigned int recover;	// number of restored symbols
} raptor_coder;


/*
 It returns a new encoder configured with given parameters.
The number of source symbols must be range of 4 ~ 8192.
The symbol size does not need to be a multiple of 8.
*/
raptor_coder *raptor_encoder_new(
	unsigned int K,		// number of source symbols
	unsigned int T);	// size of each symbol in bytes

/*
 It frees up any resources used by a decoder/encoder.
After all tasks were done, call this function.
*/
void raptor_free(raptor_coder *rc);

/*
 It adds some source symbols to the encoder.
 It returns 0 at success, see RAPTOR_ADD_* for others.
 It keeps list of added symbols.
If same source symbols are added, later one is ignored.
When it got enough symbols, it generates intermediate symbols automatically.
*/
int raptor_encoder_add(
	raptor_coder *rc,
	void *data,				// input bytes of some adding source symbols
	unsigned int off_id,	// offset to the first Encoding Symbol ID
	unsigned int num_sym);	// number of adding source symbols

/*
 It creates some encoding symbols at once.
When ESI < number of source symbols, source symbols are restored.
When ESI >= number of source symbols, repair symbols are created.
 It returns 0 at success, see RAPTOR_ENC_* for others.
*/
int raptor_encoder_create(
	raptor_coder *rc,
	void *data,				// output bytes of some created symbols
	unsigned int off_id,	// offset to the first encoding symbol ID
	unsigned int num_sym);	// number of creating symbols

/*
 It returns a new decoder initialized with given parameters.
These parameters must be same as used values in encoder.
 max_esi is the max ESI of possible repair symbols (should be larger than K).
Set the largest Encoding Symbol ID allowed at decoding.
If number of repair symbols is R, their ESIs are K ~ K+R-1.
When total number of source symbols and repair symbols is K+R = N, max_esi is K+R-1 = N-1.
*/
raptor_coder *raptor_decoder_new(
	unsigned int K,		// number of source symbols
	unsigned int T,		// size of each symbol in bytes
	unsigned int max_esi);

/*
 It returns how many extra symbols may be required to restore whole symbols.
This isn't a determined value, but is just an estimated value.
Adding more extra symbols than this value would be faster to decode.
*/
unsigned int raptor_decoder_overhead(raptor_coder *rc);

/*
 It adds an encoding symbol to the decoder.
 It returns 0 at success, see RAPTOR_ADD_* for others.
 It distinguishes the added encoding symbol by its Encoding Symbol ID.
When adding a source symbol, its ESI is range of 0 ~ K-1.
When adding a repair symbol, its ESI is range of K ~ 65535.
*/
int raptor_decoder_add(
	raptor_coder *rc,
	void *data,				// input bytes of an adding symbol
	unsigned int esi);		// encoding symbol ID

/*
 It restores intermediate symbols from added symbols in decoder.
 It returns 0 at success, see RAPTOR_ENC_* for others.
When added symbols were not enough, try again after adding more symbols.
 This tries to solve by Peeling Algorithm.
It's very fast (linear time), but would require many (several hundred) overheads.
*/
int raptor_decoder_solve_iterative(raptor_coder *rc);

/*
 It restores intermediate symbols from added symbols in decoder.
 It returns 0 at success, see RAPTOR_ENC_* for others.
When added symbols were not enough, try again after adding more symbols.
 This tries to solve by Hybrid Decoding (Peeling and Gaussian Elimination).
Though it's not so fast, it would accomplish with a few (upto ten) overheads.
*/
int raptor_decoder_solve_hybrid(raptor_coder *rc);

/*
 It restores some encoding symbols in the decoder.
When ESI < number of source symbols, source symbols are recovered.
When ESI >= number of source symbols, repair symbols are restored.
 It returns 0 at success, see RAPTOR_ENC_* for others.
*/
int raptor_decoder_recover(
	raptor_coder *rc,
	void *data,				// output bytes of some recovered symbols
	unsigned int off_id,	// offset to the first encoding symbol ID
	unsigned int num_sym);	// number of recovering symbols

#endif
