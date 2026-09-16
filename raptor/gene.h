#ifndef RAPTOR_GENERATOR_H
#define RAPTOR_GENERATOR_H

unsigned int Deg(unsigned int v);

void Trip(
	unsigned int K,
	unsigned int Ld,
	unsigned int X,		// input X = Encoding Symbol ID
	unsigned int *d, unsigned int *a, unsigned int *b);	// output (d, a, b)

#endif
