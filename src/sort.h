#ifndef SORT_HEADER
#define SORT_HEADER

unsigned int sort(FILE* restrict file, FILE* restrict goodreads, FILE* restrict candidatereads, FILE* restrict badreads,
          char* restrict head, char* restrict tail, unsigned int headlength, unsigned int taillength,
          unsigned int offset, unsigned int librarylength, unsigned int gooderrors, unsigned int maxerrors,
          unsigned int* restrict headlocationstatistics, unsigned int* restrict taillocationstatistics,
          bool secondary, bool revlibrary);

#endif
