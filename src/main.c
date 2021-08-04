#include"main.h"
#include"structures.h"
#include"sort.h"
#include"hashmap.h"
#include"print.h"
#include"mergesort.h"
#include"locationstatistics.h"

unsigned int prime = 1782629; //arbitrary large prime number. I recommend something within an order of magnitude of the (largest) sample size you intend you run.
hashmap *good, *candidate;
FILE *runinfo;
bool allnormalkeys = true; //this needs to be bugtested (it SHOuLD work...)

int main(int argc, char *argv[]){

  //overall clock start
  clock_t overallstart, overallend;
  overallstart = clock();

  //some variables that you may want to modify and recompile
  unsigned int errorallowance = 0; //errors allowed in secondary candidate sorts (+1/-1 lib length), default 0


  //mini user-manual/help
  if (argc == 1) {
    printf("Sorting using default parameters.\n");
    printf("Default parameters: sorted data.txt 0 15 ACAC CACA 0 1 false\n");
  } else if ((argc == 2) && ((strcmp(argv[1], "help") == 0) || (strcmp(argv[1], "h") == 0))){
    printf("This script sorts DNA sequence reads into good, candidate, and bad reads.\n");
    printf("This script assumes that you have preprocessed the sequencing data such that each line only contains 1 sequence.\n");
    printf("You may run the preprocess script to perform preprocessing with: ./preprocess\n");
    printf("Please have parameters following the script in the following format:\n");
    printf("./sort jobname datafile offset librarylength head tail gooderrors maxerrors revlibrary\n");
    printf("You can choose to run this script with the default parameters with:\n");
    printf("./sort\n");
    printf("Default: sorted data.txt 0 15 ACAC CACA 0 1 false\n");
    return 0;
  } else if (argc == 10) {
    printf("Running job: \"%s\" on \"%s\", with offset: %s, librarylength: %s, head: \"%s\", tail: \"%s\", gooderrors: %s, maxerrors: %s, reverse library: %s.\n", argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9]);
  } else {
    printf("There appears to be something wrong with your input parameters.\n");
    printf("Please have parameters following the script in the following format:\n");
    printf("./sort jobname datafile offset librarylength head tail gooderrors maxerrors revlibrary.\n");
    return 0;
  }

  //declare some variables
  char *job, *inputfile, *head, *tail, *librarylengthstring;
  unsigned int offset, librarylength, gooderrors, maxerrors;
  unsigned int joblen;
  char *outputdirectory, *goodreadsfp, *goodreadscountfp, *candidatereadsfpbase, *candidatereadsfp, *candidatereadscountfp, *badreadsfp, *totalcountfp;
  size_t outputdirectorylen, goodreadsfplen, goodreadscountfplen, candidatereadsfpbaselen, candidatereadsfplen, candidatereadscountfplen, badreadsfplen, totalcountfplen;
  bool revlibrary;

  //initialize some variables to either default or user-input
  if (argc == 1){ //default
    job = "sorted";
    inputfile ="data.txt";
    head = "ACAC";
    tail = "CACA";
    offset = 0;
    librarylengthstring = "15";
    librarylength = 15;
    gooderrors = 0;
    maxerrors = 1;
    revlibrary = false;
  } else if (argc == 10){  //user-input
    job = argv[1];
    inputfile = argv[2];
    offset = (unsigned int) strtol(argv[3], NULL, 10);
    librarylengthstring = argv[4];
    librarylength = (unsigned int) strtol(librarylengthstring, NULL, 10);
    head = argv[5];
    tail = argv[6];
    gooderrors = (unsigned int) strtol(argv[7], NULL, 10);
    maxerrors = (unsigned int) strtol(argv[8], NULL, 10);
    if (gooderrors >= maxerrors){
      printf("Good errors must be less than max errors.\n");
      return 0;
    }
    if (strcmp(argv[9], "true") == 0){
      revlibrary = true;
    } else if (strcmp(argv[9], "false") == 0){
      revlibrary = false;
    } else {
      printf("Incorrect parameter for revlibrary.\n");
      return 0;
    }
  } else { //this should not be reached
    printf("Something went horribly wrong since you should never have reached this line. Type \"./sort help\" for help and contact the developer.\n");
    return 0;
  }

  //some dynamic strings for filepaths
  joblen = strlen(job);

  outputdirectorylen = joblen + 3; //+2 for "./", +1 for '\0';
  outputdirectory = malloc(outputdirectorylen * sizeof(char));
  snprintf(outputdirectory, outputdirectorylen, "./%s", job);
  mkdir(outputdirectory, 0777); //does not error check

  goodreadsfplen = joblen + joblen + 16; //+13 for "goodreads.txt", +1 for '/', +1 for '_', +1 for '\0'
  goodreadsfp = malloc(goodreadsfplen * sizeof(char));
  snprintf(goodreadsfp, goodreadsfplen, "%s/%s_goodreads.txt", job, job);

  goodreadscountfplen = joblen + joblen + 21; //+18 for "goodreadscount.txt", +1 for '/', +1 for '_', +1 for '\0'
  goodreadscountfp = malloc(goodreadscountfplen * sizeof(char));
  snprintf(goodreadscountfp, goodreadscountfplen, "%s/%s_goodreadscount.txt", job, job);

  candidatereadsfpbaselen = joblen + joblen + 18; //+14 for "candidatereads", +1 for '/', +2 for '_', +1 for '\0'
  candidatereadsfpbase = malloc(candidatereadsfpbaselen * sizeof(char));
  snprintf(candidatereadsfpbase, candidatereadsfpbaselen, "%s/%s_candidatereads_", job, job);

  candidatereadsfplen = candidatereadsfpbaselen + strlen(librarylengthstring) + 4; //+4 for ".txt", already includes +1 for '\0'
  candidatereadsfp = malloc(candidatereadsfplen * sizeof(char));
  snprintf(candidatereadsfp, candidatereadsfplen, "%s%s.txt", candidatereadsfpbase, librarylengthstring);

  candidatereadscountfplen = joblen + joblen + 26; //+23 for "candidatereadscount.txt", +1 for '/', +1 for '_', +1 for '\0'
  candidatereadscountfp = malloc(candidatereadscountfplen * sizeof(char));
  snprintf(candidatereadscountfp, candidatereadscountfplen, "%s/%s_candidatereadscount.txt", job, job);

  badreadsfplen = joblen + joblen + 15; //+12 for "badreads.txt", +1 for '/', +1 for '_', +1 for '\0'
  badreadsfp = malloc(badreadsfplen * sizeof(char));
  snprintf(badreadsfp, badreadsfplen, "%s/%s_badreads.txt", job, job);

  totalcountfplen = joblen + joblen + 17; //+14 for "totalcount.txt", +1 for '/', +1 for '_', +1 for '\0'
  totalcountfp = malloc(totalcountfplen * sizeof(char));
  snprintf(totalcountfp, totalcountfplen, "%s/%s_totalcount.txt", job, job);

  //script variables
  unsigned int i;
  unsigned int headlength, taillength;
  headlength = strlen(head);
  taillength = strlen(tail);

  //for locationstatistics
  unsigned int headlocationstatistics[headlength * 5];
  unsigned int taillocationstatistics[taillength * 5];
  for(i = 0; i < headlength * 5; i++){
    headlocationstatistics[i] = 0;
  }
  for(i = 0; i < taillength * 5; i++){
    taillocationstatistics[i] = 0;
  }

  FILE *file = fopen(inputfile, "r");
  if(!file){
    printf("Error opening input data file.\n");
    return 0;
  }

  FILE *badreads = fopen(badreadsfp , "w");
  if(!badreads){
    printf("Error creating bad reads file for write.\n");
    return 0;
  }

  FILE *goodreads = fopen(goodreadsfp, "w");
  if(!goodreads){
    printf("Error opening %s for write.\n", goodreadsfp);
    return 0;
  }

  FILE *candidatereads = fopen(candidatereadsfp, "w");
  if(!candidatereads){
    printf("Error opening %s for write.\n", candidatereadsfp);
    return 0;
  }

  //some dynamic memory allocation
  good = (hashmap*) malloc(prime * sizeof(hashmap*));
  candidate = (hashmap*) malloc(prime * sizeof(hashmap*));

  //initialize hashmaps w/ linked lists
  initialize(good);
  initialize(candidate);

  //begin initial sorting
  printf("Opening \"%s\" for reading...\n",inputfile);

  clock_t start = clock();

  //close files
  sort(file, goodreads, candidatereads, badreads, head, tail, headlength, taillength,
            offset, librarylength, gooderrors, maxerrors,
            headlocationstatistics, taillocationstatistics, false, revlibrary);
  fclose(file);
  fclose(badreads);

  clock_t end = clock();

  printf("Finished initial sort of \"%s\" in %f seconds.\n",inputfile, (double)(end - start) / CLOCKS_PER_SEC);

  //printlocationstatistics
  printlocationstatistics(outputdirectory, headlocationstatistics, taillocationstatistics, head, tail);

  //free output directory as it is no longer needed
  free(outputdirectory);

  //runinfo
  char RIfile[200];
  snprintf(RIfile, 199, "./%s/runinfo.txt", job);

  runinfo = fopen(RIfile, "w");
  if (!runinfo) {
    printf("Error opening %s for write.\n", RIfile);
    return 0;
  }

  fprintf(runinfo, "Running job: \"%s\" on \"%s\", with offset: %u, librarylength: %u, head: \"%s\", tail: \"%s\", gooderrors: %u, maxerrors: %u.\n\n", job, inputfile, offset, librarylength, head, tail, gooderrors, maxerrors);

  printf("Processing initial sort...\n");
  start = clock();

  char tempfile[200];
  snprintf(tempfile, 199, "%s/%s_badreadstemp.txt", job, job);
  FILE *badreadstemp = fopen(tempfile, "w");

  badreads = fopen(badreadsfp, "r");

  if (offset > 0) {
    //note that since we use gooderrors for maxerrors as well, this should never send anything to candidate
    sort(badreads, goodreads, candidatereads, badreadstemp, head, tail, headlength, taillength,
            offset-1, librarylength, gooderrors, gooderrors,
            headlocationstatistics, taillocationstatistics, false, revlibrary);
  } else if (offset == 0){
    char newhead[headlength];
    for(i = 0; i<headlength-1; i++){
      newhead[i] = head[i+1];
    }
    newhead[headlength] = '\0'; //not sure if this is needed
    //note with truncated head, the locationstatistics will be wrong if gooderrors > 0
    sort(badreads, goodreads, candidatereads, badreadstemp, newhead, tail, headlength-1, taillength,
            offset, librarylength, gooderrors, gooderrors,
            headlocationstatistics, taillocationstatistics, false, revlibrary);
  } else {
    printf("How did an unsigned int not end up being greater than or equal to 0?\n");
  }

  fclose(badreads);
  fclose(candidatereads);
  fclose(badreadstemp);
  remove(badreadsfp);
  rename(tempfile, badreadsfp);

  //connect the linked lists in good hashmap
  node *temp, *prev;
  node *goodheadnode;
  goodheadnode = linkentirelist(good);

  splitlist twolists = split_list(goodheadnode); //split the list for mergesort since we expect most values (>>95%) to be n=1
  goodheadnode = twolists.list2; //mergesort only on the values that are not n=1
  mergesort(&goodheadnode); //iterative, result is in descending order now
  temp = goodheadnode;

  //this is necessary if list2 (n>1 counts) is NULL, since append order is list2 -> list 1. Obviously doesn't matter if list 1 is null.
  while(temp){
    prev=temp;
    temp=temp->next;
  }
  if(goodheadnode){
    prev->next = twolists.list1;
  } else {
    goodheadnode = twolists.list1;
  }

  initialize(good); //necessary since the good[index].head nodes still point to certain nodes //could probably remove that pointer during linkentirelist
  //print goodreads count and free linked list and free good
  unsigned int goodcount = printgoodcount(goodheadnode, goodreadscountfp, runinfo);

  end = clock();
  printf("Processing of initial sort finished in %f seconds.\n", (double)(end - start) / CLOCKS_PER_SEC);

  printf("Starting secondary sort...\n");
  //catch candidates that are +1/-1 of library length
  unsigned int newlibrarylength = librarylength - 1;
  char newlibrarylengthbuffer[4];
  char *newlibrarylengthstring;
  snprintf(newlibrarylengthbuffer, 3, "%u", newlibrarylength);
  newlibrarylengthstring = newlibrarylengthbuffer;
  if (strlen(newlibrarylengthstring) != strlen(librarylengthstring)) {
    candidatereadsfplen = candidatereadsfpbaselen + strlen(newlibrarylengthstring) + 4; //+4 for ".txt"
    candidatereadsfp = realloc(candidatereadsfp, candidatereadsfplen);
  }
  snprintf(candidatereadsfp, candidatereadsfplen, "%s%s.txt", candidatereadsfpbase, newlibrarylengthstring);
  candidatereads = fopen(candidatereadsfp, "w");
  if (!candidatereads){
    printf("Error opening %s for write.\n", candidatereadsfp);
    return 0;
  }
  badreads = fopen(badreadsfp, "r");
  badreadstemp = fopen(tempfile, "w");
  sort(badreads, goodreads, candidatereads, badreadstemp, head, tail, headlength, taillength,
            offset, librarylength-1, gooderrors, maxerrors+errorallowance,
            headlocationstatistics, taillocationstatistics, true, revlibrary);

  fclose(candidatereads);
  fclose(badreads);
  fclose(badreadstemp);

  printf("Continuing secondary sort...\n");

  //candidate catcher part 2
  newlibrarylength = librarylength + 1;
  snprintf(newlibrarylengthbuffer, 3, "%d", newlibrarylength); //again new librarylength as string
  newlibrarylengthstring = newlibrarylengthbuffer;
  if ((candidatereadsfpbaselen + strlen(newlibrarylengthstring) + 4) != candidatereadsfplen) {
    candidatereadsfplen = candidatereadsfpbaselen + strlen(newlibrarylengthstring) + 4; //+4 for ".txt"
    candidatereadsfp = realloc(candidatereadsfp, candidatereadsfplen);
  }
  snprintf(candidatereadsfp, candidatereadsfplen, "%s%s.txt", candidatereadsfpbase, newlibrarylengthstring);
  candidatereads = fopen(candidatereadsfp, "w");
  badreads = fopen(badreadsfp, "w");
  badreadstemp = fopen(tempfile, "r");

  unsigned int badcount = sort(badreadstemp, goodreads, candidatereads, badreads, head, tail, headlength, taillength,
            offset, librarylength+1, gooderrors, maxerrors+errorallowance,
            headlocationstatistics, taillocationstatistics, true, revlibrary);

  fclose(goodreads); //close goodreads now that we're done sorting.
  fclose(candidatereads);
  fclose(badreads);
  fclose(badreadstemp);
  remove(tempfile);

  printf("Secondary sort finished.\n");

  //print candidate reads count and add candidate nodes into good for a total
  FILE *candidatereadscount = fopen(candidatereadscountfp, "w");
  node *candidateheadnode;
  candidateheadnode = linkentirelist(candidate);
  printf("Processing secondary sort...\n");
  mergesort(&candidateheadnode);
  temp = candidateheadnode;
  unsigned int candidatecount = 0;
  while (temp != NULL){
    candidatecount = candidatecount + temp->count;
    fprintf(candidatereadscount, "%s\t%u\n", temp->library,temp->count);
    prev=temp;
    temp = move(good, temp); //move returns next node in list, and NULL if no nodes remain
  }

  fclose(candidatereadscount);
  free(candidate); //free candidate now that we're done with it

  printf("Processing total count...\n");
  start = clock();

  node *totalheadnode;
  totalheadnode = linkentirelist(good);
  twolists = split_list(totalheadnode);
  totalheadnode = twolists.list2;

  mergesort(&totalheadnode);
  temp = totalheadnode;

  //this is necessary if list2 (the >1 counts) is NULL since append order is list2 then list1
  while(temp){
    prev = temp;
    temp = temp->next;
  }
  if(totalheadnode){
    prev->next = twolists.list1;
    temp=totalheadnode;
  } else {
    temp = twolists.list1;
  }

  //iterate through linkedlist and print to file
  FILE *totalcountfile = fopen(totalcountfp, "w");
  while (temp){
    fprintf(totalcountfile, "%s\t%u\n", temp->library, temp->count);
    prev = temp;
    temp = temp->next;
    //free nodes in good now that we're done with them
    free(prev->library);
    free(prev);
  }

  free(good); //free good now that we're done with it
  fclose(totalcountfile);

  end = clock();
  printf("Total count processing took: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);

  //print some stuff to runinfo
  fprintf(runinfo, "\nTotal Number of Reads: %u\n", goodcount+candidatecount+badcount);
  fprintf(runinfo, "Number of Good Reads: %u\n", goodcount);
  fprintf(runinfo, "Number of Candidate Reads: %u\n", candidatecount);
  fprintf(runinfo, "Number of Bad Reads: %u\n", badcount);

  //close runinfo
  fclose(runinfo);

  unsigned int totalcount;
  float goodpercent = 0.00, candidatepercent = 0.00, badpercent = 0.00;
  totalcount = goodcount + candidatecount + badcount;
  FILE *exportfile;
  char *exportfp;
  size_t exportfplen;
  exportfplen = strlen(job) + 12; //+1 for '\0', +1 for '/', +10 for "export.csv"
  exportfp = malloc(exportfplen * sizeof(char));
  snprintf(exportfp, exportfplen, "%s/export.csv", job);
  exportfile = fopen(exportfp, "w");
  if (!exportfile) {
    printf("Error opening %s for write.\n", exportfp);
    return 0;
  }

  if (totalcount != 0) {
    goodpercent = 100 * (float) goodcount / (float) totalcount;
    candidatepercent = 100 * (float) candidatecount / (float) totalcount;
    badpercent = 100 * (float) badcount / (float) totalcount;
  }

  fprintf(exportfile, "%s\t%u\t%u\t%u\t%u\t%.2f\t%.2f\t%.2f\t%u\t%u\t%s\t%s\t%u\t%u", job, totalcount, goodcount, candidatecount, badcount, goodpercent, candidatepercent, badpercent, librarylength, offset, head, tail, gooderrors, maxerrors);
  fclose(exportfile);

  //free
  free(goodreadsfp);
  free(goodreadscountfp);
  free(candidatereadsfpbase);
  free(candidatereadsfp);
  free(candidatereadscountfp);
  free(badreadsfp);
  free(totalcountfp);
  free(exportfp);

  //end overall clock
  overallend = clock();

  printf("Sorting finished in: %f seconds. Please check the results in the \"%s\" directory.\n", (double)(overallend - overallstart) / CLOCKS_PER_SEC, job);

  //keep in mind that since we rely on the key to get the index, if we have both minikeys and longkeys, the index they are inserted to may vary.

  return 1;
}
