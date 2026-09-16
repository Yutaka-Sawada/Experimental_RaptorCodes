# Experimental_RaptorCodes

This is an experimental implementation of Raptor Codes. 
It only supports core feature (encoder and decoder) of RFC 5053.

The implementation consists of files in `raptor` directory. 
It's possible to setup encoder/decoder, some symbols, and test recovery. 
For the usage, refer `raptor.h`. 
I put a sample `main.c` to use the Recovery Codes.

## Multiple decoding methods

For recovery, I implemented Peeling Algorithm, Gaussian Elimination, and Inactivation Decoding.
- Peeling Algorithm (called as Peeling Decoder or Iterative Decoding) is the fastest, but requires many overheads.
- Gaussian Elimination is slow, but requires less overheads.
- Inactivation Decoding is faster than Gaussian Elimination. 
Because Inactivation Decoding has patent issue, it's disabled by default. 
If you want to test their speed, you need to edit source code.

## Caution

Because I didn't compare result of encoded symbols, 
there may be a compatibility issue with other implementations of Raptor Codes.

If you want to use this for practical usage, you need to check and test by yourself.

