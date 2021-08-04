#include"main.h"
#include"structures.h"
#include"print.h"
#include"hashmap.h"

extern hashmap *good;
//extern FILE *runinfo;

void printlocationstatistics(char* restrict outputdirectory, unsigned int* restrict headlocationstatistics, unsigned int* restrict taillocationstatistics, char* restrict head, char* restrict tail){ //print locationstatistics

  //declare some variables
  //note that the performlocatsionstatistics awkwardly made i = row and j = col
  unsigned int i, j;
  unsigned int headlength = strlen(head);
  unsigned int taillength = strlen(tail);
  char *nucleotides = "ATCGN"; //this can be changed

  //get filepath and open
  char *locationstatisticsfp;
  size_t locationstatisticsfplen;
  FILE *locationstatistics;
  locationstatisticsfplen = strlen(outputdirectory) + 24; //+22 for "locationstatistics.txt", +1 for '/', +1 for '\0'
  locationstatisticsfp = malloc(locationstatisticsfplen * sizeof(char));
  snprintf(locationstatisticsfp, locationstatisticsfplen, "%s/locationstatistics.txt", outputdirectory);
  locationstatistics = fopen(locationstatisticsfp, "w");

  if (!locationstatistics) {
    printf("Error opening %s for write.\n", locationstatisticsfp);
    return;
  }

  //free filepath string
  free(locationstatisticsfp);

  //print head nucleotides as column name
  fputs("*",locationstatistics);
  for(i = 0; i < headlength; i++){
    fprintf(locationstatistics, "\t%c", head[i]);
  }
  fputs("\n",locationstatistics);

  //print the matrix for the head
  for(i = 0; i < 5; i++){
    fprintf(locationstatistics,"%c",nucleotides[i]);
    for(j = 0; j < headlength; j++){
      fprintf(locationstatistics,"\t%d",headlocationstatistics[i + j*5]);
    }
    fputs("\n", locationstatistics);
  }
  fputs("\n*",locationstatistics);

  //print tail nucleotides as column name
  for(i = 0; i < taillength; i++){
    fprintf(locationstatistics, "\t%c", tail[i]);
  }
  fputs("\n",locationstatistics);

  //print the matrix for the tail
  for(i = 0; i < 5; i++){
    fprintf(locationstatistics,"%c",nucleotides[i]);
    for(j = 0; j < taillength; j++){
      fprintf(locationstatistics,"\t%d",taillocationstatistics[i + j*5]);
    }
    fputs("\n", locationstatistics);
  }

  //close file
  fclose(locationstatistics);

  return;
}

unsigned int printgoodcount(node* goodheadnode, char* restrict goodreadscountfp, FILE* restrict runinfo){
  if (!goodheadnode){
    printf("Good head node is NULL.\n");
    return 0;
  }

  FILE *goodreadscount = fopen(goodreadscountfp,"w");
  if (!goodreadscount){
    printf("Error opening %s for write.\n", goodreadscountfp);
    return 0;
  }

  node *temp;
  temp = goodheadnode;

  unsigned int goodcount = 0; //total good reads
  unsigned int prevcount = temp->count;

  unsigned int goodtotal = 0; //number of good reads for count=n
  unsigned int goodnodes = 0; //total number of good read nodes
  unsigned int lastcount;

  listnode *listhead = (listnode*) malloc(sizeof( listnode));
  listnode *listtemp = NULL;

  //write to file and create new linked list for count percents while freeing nodes for memory
  while (temp != NULL) { //continues until end of list

    //linked list for count percents
    if (prevcount == temp->count){ //if library sequence has same count as previous
      goodtotal++;
    } else { //if current library sequence has a different count
      if (listtemp){ //append new node to linked list (for previous count)
        listnode *listnew = ( listnode*) malloc(sizeof( listnode));
        listnew->count = prevcount;
        listnew->totalcount = goodtotal;
        listnew->next = NULL;
        listtemp->next = listnew;
        listtemp = listnew;
      } else { //create head node (for previous count)
        listhead->count = prevcount;
        listhead->totalcount = goodtotal;
        listhead->next = NULL;
        listtemp = listhead;
      }
      goodnodes = goodnodes + goodtotal;
      goodtotal = 1; //set good total to 1 for current count
    }
    prevcount = temp->count; //update prevcount to current sequence's count

    //get total number of good reads for runinfo
    goodcount = goodcount + temp->count;

    //write to good reads count file
    //we should change this to be a tab delimited file (col1 = sequence, col2 = count)
    fprintf(goodreadscount, "%s\t%u\n", temp->library, temp->count);
    //fprintf(goodreadscount,"%s\tn=%d\n", temp->library, temp->count);

    //iterate down list while freeing
    lastcount = temp->count;
    temp = move(good, temp); //moves nodes back into the good hashmap for a total count with candidates
  }
  goodnodes = goodnodes + goodtotal;

  //this part is necessary since the above loop runs for the previous count (and thus the current [last] count needs to be added)
  if (listtemp){
     listnode *listnew = ( listnode*) malloc(sizeof( listnode));
    listnew->count = lastcount;
    listnew->totalcount = goodtotal;
    listnew->next = NULL;
    listtemp->next = listnew;
    listtemp = listnew;
  } else {
    listhead->count = lastcount;
    listhead->totalcount = goodtotal;
    listhead->next = NULL;
    listtemp = listhead;
  }

  //close file
  fclose(goodreadscount);

  //print the count percents to runinfo and free linked list
  fputs("Good Reads Count Values given as a total and a percentage of total good reads:\n", runinfo);

  //declare some variables
  double percent;
  listtemp = listhead;

  //prints to runinfo then frees current node
  while (listtemp){
    listhead = listtemp;
    listtemp = listtemp->next;
    percent = 100 * (double) listhead->totalcount / (double) goodnodes;
    fprintf(runinfo, "Count: %d\tTotal: %d\tPercentage: %.5f%%\n", listhead->count, listhead->totalcount, percent);
    free(listhead);
  }

  return goodcount;
}
