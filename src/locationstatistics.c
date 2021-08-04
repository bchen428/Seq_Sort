#include"main.h"
#include"locationstatistics.h"

//updates locationstatistics matrices
void performlocationstatistics(unsigned int* restrict headlocationstatistics, unsigned int* restrict taillocationstatistics,
  char* restrict head, char* restrict tail, char* restrict headseq, char* restrict tailseq,
  unsigned int headlength, unsigned int taillength){

  //declare some variables
  unsigned int i;
  unsigned int numheadseq[headlength], numtailseq[taillength];

  //loop for head
  for(i = 0; i < headlength; i++){ //convert library sequence to numeric
    if(headseq[i] == 'A'){
      numheadseq[i] = '0';
    } else if (headseq[i] == 'T'){
      numheadseq[i] = '1';
    } else if (headseq[i] == 'C'){
      numheadseq[i] = '2';
    } else if (headseq[i] == 'G'){
      numheadseq[i] = '3';
    } else { //currently sets all non ATCG to 4.
      numheadseq[i] = '4';
    }
    if(head[i] != headseq[i]){
      headlocationstatistics[i*5 + (numheadseq[i]-'0')]++;
    }
  }

  //loop for tail
  for(i = 0; i < taillength; i++){ //convert library sequence to numeric
    if(tailseq[i] == 'A'){
      numtailseq[i] = '0';
    } else if (tailseq[i] == 'T'){
      numtailseq[i] = '1';
    } else if (tailseq[i] == 'C'){
      numtailseq[i] = '2';
    } else if (tailseq[i] == 'G'){
      numtailseq[i] = '3';
    } else { //currently sets all non ATCG to 4.
      numtailseq[i] = '4';
    }
    if(tail[i] != tailseq[i]){
      taillocationstatistics[i*5 + (numtailseq[i]-'0')]++;
    }
  }
}
