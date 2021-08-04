#include"main.h"
#include"sort.h"
#include"structures.h"
#include"hashmap.h"
#include"locationstatistics.h"

extern hashmap *good, *candidate;

unsigned int sort(FILE* restrict file, FILE* restrict goodreads, FILE* restrict candidatereads, FILE* restrict badreads,
          char* restrict head, char* restrict tail, unsigned int headlength, unsigned int taillength,
          unsigned int offset, unsigned int librarylength, unsigned int gooderrors, unsigned int maxerrors,
          unsigned int* restrict headlocationstatistics, unsigned int* restrict taillocationstatistics,
          bool secondary, bool revlibrary){

  //initial sort that takes the preprocessed fastq data and gives goodreads, candidatereads, and badreads
  size_t line_size = 200; // chose 200 arbitrarily. can increase this if necessary. Yes, we could opt to read the first line and count the characters... but... 200 should be more than enough.
  size_t library_size;
  char *sequence = (char*) malloc(line_size * sizeof(char));

  unsigned int badcount = 0;

  unsigned int headendindex = offset + headlength;
  unsigned int tailindex = librarylength + headendindex;
  unsigned int i, j;
  char headseq[headlength], tailseq[taillength];
  char restofseq[(200-tailindex-taillength)];
  char *library;
  bool headenderror, tailenderror;
  unsigned int errors, headenderrorcount, tailenderrorcount;
  u_int64_t key;
  unsigned int remaining;
  unsigned int notzero = offset;
  if(offset == 0){
    notzero = 1;
  }
  char offsetseq[notzero];
  offsetseq[0] = '\0';
  char c;

  while(getline(&sequence, &line_size, file) != -1){ //continues as long as there are lines in file

    if(secondary){ //this part should be changed to suit needs
      errors = 1;
    } else {
      errors = 0;
    }

    headenderror = true;
    tailenderror = true;
    headenderrorcount = 0;
    tailenderrorcount = 0;

    for(i = headlength; i > 0; i--){ //determines number of mismatches between head and sequence adapter
      if(sequence[i-1+offset] != head[i-1]){ //(starts from right side of head so as to count possible deletions)
        errors++;
        if(headenderror == true){
          headenderrorcount++;
        }
      } else {
        headenderror = false; //end of the contiguous errors
      }
    }

    for (i = 0; i<taillength; i++){ //determines number of mismatches between tail and sequence adapter
      if(sequence[i+tailindex] != tail[i]){ //starts from left side
        errors++;
        if(tailenderror == true){
          tailenderrorcount++;
        }
      } else {
        tailenderror = false; //end of the continguous errors
      }
    }

    library_size = librarylength + headenderrorcount + tailenderrorcount + 1; //+1 for '\0'
    library = (char*) malloc(library_size * sizeof(char));

    if(secondary){
      if((headenderrorcount > 0) != (tailenderrorcount > 0)){
        errors--;
      }
    }

    if (errors <= gooderrors){ //inserts to good   [unsigned so cannot be negative]
      //for locationstatistics
      strncpy(headseq, sequence+offset, headlength);
      headseq[headlength]='\0';
      strncpy(tailseq, sequence+tailindex, taillength);
      tailseq[taillength]='\0';
      performlocationstatistics(headlocationstatistics, taillocationstatistics, head, tail, headseq, tailseq, headlength, taillength); //update locationstatistics matrix

      strncpy(offsetseq, sequence, offset);
      offsetseq[offset]='\0';
      remaining = strlen(sequence)-tailindex-taillength;
      strncpy(restofseq, sequence+tailindex+taillength, remaining);
      restofseq[remaining]='\0';
      if (headenderrorcount > 0 && tailenderrorcount > 0){ //if deletion from both ends
        strncpy(library, sequence+headendindex-headenderrorcount, librarylength+headenderrorcount+tailenderrorcount);
        library[librarylength+headenderrorcount+tailenderrorcount]='\0'; //this may require malloc?
        for(i = 1; i<=headenderrorcount; i++){
          headseq[headlength-i] = '*';
        }
        for(i = 0; i<tailenderrorcount; i++){
          tailseq[i]='*';
        }
      } else if (headenderrorcount > 0){ //if deletion from head
        strncpy(library, sequence+headendindex-headenderrorcount, librarylength+headenderrorcount); //note uses strncpy
        library[librarylength+headenderrorcount]='\0';
        for(i = 1; i<=headenderrorcount; i++){
          headseq[headlength-i] = '*';
        }
      } else if (tailenderrorcount > 0){ //if deletion from tail
        strncpy(library, sequence+headendindex, librarylength+tailenderrorcount); //note uses strncpy
        library[librarylength+tailenderrorcount]='\0';
        for(i = 0; i<tailenderrorcount; i++){
          tailseq[i]='*';
        }
      } else { //no deletions
        strncpy(library, sequence+headendindex, librarylength); //note uses strncpy
        library[librarylength]='\0';
      }
      if (revlibrary){
        j = librarylength/2;
        for(i = 0; i < j; i++){
          c = library[i];
          library[i] = library[librarylength - i - 1];
          library[librarylength - i - 1] = c;
        }
      }
      fprintf(goodreads, "%s\t%s\t%s\t%s\t%s", offsetseq, headseq, library, tailseq, restofseq);
      key = library_to_key(library, librarylength+headenderrorcount+tailenderrorcount);
      insert(good, key, library);

    } else if (errors <= maxerrors){ //inserts to candidate [checks the first condition so should]

      //for locationstatistics
      strncpy(headseq, sequence+offset, headlength);
      headseq[headlength]='\0';
      strncpy(tailseq, sequence+tailindex, taillength);
      tailseq[taillength]='\0';
      performlocationstatistics(headlocationstatistics, taillocationstatistics, head, tail, headseq, tailseq, headlength, taillength); //update locationstatistics matrix

      strncpy(offsetseq, sequence, offset);
      offsetseq[offset]='\0';
      remaining = strlen(sequence)-tailindex-taillength;
      strncpy(restofseq, sequence+tailindex+taillength, remaining);
      restofseq[remaining]='\0';
      if(secondary){
        if (headenderrorcount > 0 && tailenderrorcount > 0){ //if deletion from both ends
          strncpy(library, sequence+headendindex-headenderrorcount, librarylength+headenderrorcount+tailenderrorcount);
          library[librarylength+headenderrorcount+tailenderrorcount]='\0';
          for(i = 1; i<=headenderrorcount; i++){
            headseq[headlength-i] = '*';
          }
          for(i = 0; i<tailenderrorcount; i++){
            tailseq[i]='*';
          }
        } else if (headenderrorcount > 0){ //if deletion from head
          strncpy(library, sequence+headendindex-headenderrorcount, librarylength+headenderrorcount); //note uses strncpy
          library[librarylength+headenderrorcount]='\0';
          for(i = 1; i<=headenderrorcount; i++){//doublecheckthis
            headseq[headlength-i] = '*';
          }
        } else if (tailenderrorcount > 0){ //if deletion from tail
          strncpy(library, sequence+headendindex, librarylength+tailenderrorcount); //note uses strncpy
          library[librarylength+tailenderrorcount]='\0';
          for(i = 0; i<tailenderrorcount; i++){
            tailseq[i]='*';
          }
        } else {
          strncpy(library, sequence+headendindex, librarylength); //note uses strncpy
          library[librarylength]='\0';
        }
      } else {
        strncpy(library, sequence+headendindex, librarylength);
        library[librarylength]='\0';
      }
      if (revlibrary){
        j = librarylength/2;
        for(i = 0; i < j; i++){
          c = library[i];
          library[i] = library[librarylength - i - 1];
          library[librarylength - i - 1] = c;
        }
      }
      fprintf(candidatereads, "%s\t%s\t%s\t%s\t%s", offsetseq, headseq, library, tailseq, restofseq);
      if(secondary){
        key = library_to_key(library, librarylength+headenderrorcount+tailenderrorcount); //gets key
      } else {
        key = library_to_key(library, librarylength);
      }
      insert(candidate, key, library);
    } else { //write to badreads

      badcount++;
      //im too lazy to reverse the library of the badreads...
      fprintf(badreads, "%s", sequence);
      free(library);
    }
  }
  free(sequence);
  return badcount;
}
