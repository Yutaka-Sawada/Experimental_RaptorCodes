#ifndef RAPTOR_PRECODE_H
#define RAPTOR_PRECODE_H

void init_precode(raptor_coder *rc);

void buildGraySequence(unsigned int *seq, unsigned int length, unsigned int b);

unsigned int init_precode_matrix(
	raptor_coder *rc,
	unsigned char *A);		// zero fill at first

int init_precode_bit_matrix(
	raptor_coder *rc,
	uint64_t *A,			// zero fill at first
	unsigned int int_col);	// number of integers for columns in matrix A

int make_generator_matrix(raptor_coder *rc);	// Gaussian Elimination
int make_inactive_matrix(raptor_coder *rc);		// Inactivation Decoding

#endif
