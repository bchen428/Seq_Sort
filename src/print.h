#ifndef PRINT_HEADER
#define PRINT_HEADER

void printlocationstatistics(char* restrict outputdirectory, unsigned int* restrict headlocationstatistics, unsigned int* restrict taillocationstatistics, char* restrict head, char* restrict tail);
unsigned int printgoodcount(node* goodheadnode, char* restrict goodreadscountfp, FILE* restrict runinfo);

#endif
