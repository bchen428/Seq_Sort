#include"main.h"
#include"structures.h"
#include"hashmap.h"

extern unsigned int prime;
extern bool allnormalkeys;

int findkey(node* restrict list, unsigned long key, char* restrict library){ //find index of node in list with matching key
  int index = 0;
  node *temp = list;
  //checking for all normal keys here should be sufficient since the insert based on index (key % prime) should not matter as different lengths should never match
  if(allnormalkeys){
    while (temp){ //continues to last node
      if (temp->key == key){ //checks if matching
        return index;
      }
    temp = temp->next;
    index++;
    }
  } else {
    while (temp){
      if (strcmp(temp->library,library) == 0){ //does strcmp since keys can match for different strings
        return index;
      }
    temp = temp->next;
    index++;
    }
  }
  return -1; //return -1 if no nodes with matching key
}

node* get_node(node* restrict list, unsigned int index){ //returns a particular node at index (find_index) of list
	unsigned int i = 0;
	node *temp = list;
	while (i < index){ //loops until temp = node at index
		temp = temp->next;
		i++;
	}
	return temp;
}

void insert(hashmap* restrict array, u_int64_t key, char* restrict library){ //insert/update node in given hashtable with specified parameters
  unsigned long index = key % prime; //calculate index in array
  node *list = (node*) array[index].head;
  if (list == NULL){ //if array head is NULL, set head to new node
    node *new = (node*) malloc(sizeof(node));
    new->key = key;
    new->library = library; //strdup(library);
    new->count = 1;
    new->next = NULL;
    array[index].head = new;
  } else { //if array head exists
    int matching_index = findkey(list, key, library); //check if exists a node in list with matching key
    if (matching_index == -1){ //if no matching node found, insert a new node at head and
      node *new = (node*) malloc(sizeof(node));
      new->key = key;
      new->library = library; //strdup(library);
      new->count = 1;
      new->next = array[index].head; //set next to current array head node
      array[index].head=new; //set array head node to new node
    } else { //update matching node count++
      node *matching = get_node(list, matching_index);
      matching->count++;
      free(library);
    }
  }
  //rehashing and stuff goes here if necessary to implement (does not seem necessary)
}

node *move(hashmap* restrict array, node* existingnode){ //move an existing node to a new hashmap
  if(array == NULL){
    printf("Error: a NULL hashmap was passed to the move function.\n");
    return NULL;
  }
  if(existingnode == NULL){
    printf("Error: a NULL node was passed to the move function.\n");
    return NULL;
  }

  node *nextnode = NULL;
  unsigned long index = existingnode->key % prime;
  node *list = (node*) array[index].head; //get first node of array

  if (list == NULL){ //if array head is NULL
    array[index].head = existingnode;
    nextnode = existingnode->next;
    existingnode->next = NULL;
  } else {
    int find_index = findkey(list, existingnode->key, existingnode->library);
    if (find_index == -1){ //no match
      nextnode = existingnode->next;
      existingnode->next = array[index].head;
      array[index].head = existingnode;
    } else { //matched
      node *matching = get_node(list, find_index);
      matching->count = matching->count + existingnode->count;
      nextnode = existingnode->next;
      free(existingnode->library);
      free(existingnode); //frees the existing node since it is no longer needed
    }
  }
  return nextnode; //returns the next node
}

void initialize(hashmap* restrict array){ //initializes the hashtable array by setting all head nodes to NULL
  for (unsigned int i = 0; i < prime; i++){
		array[i].head = NULL;
	}
}

u_int64_t library_to_key(char* restrict library, unsigned int librarylength){ //converts library sequence to unsigned long integer (max 18 digits [since strtol is signed long and that cannot be 99[restof19digitnumber]])
  unsigned int i, j;
  u_int64_t numerickey;
  char key[librarylength+1];

  if (librarylength <= 19){ //if it will fit normally into a long integer
    for(i = 0; i < librarylength; i++){ //convert library sequence to numeric
      if(library[i] == 'A'){
        key[i] = '1';
      } else if (library[i] == 'T'){
        key[i] = '2';
      } else if (library[i] == 'C'){
        key[i] = '3';
      } else if (library[i] == 'G'){
        key[i] = '4';
      } else { //currently sets all non ATCG to 5 (This also is where any other random errors would go)
        key[i] = '5';
        printf("Non-standard nucleotide detected in library.\n");
      }
    }
    key[librarylength]='\0';
    numerickey = strtoull(key, NULL, 10); //strtol(string, &randompointer, numericbase) returns string converted to long integer
  } else { //if it needs to be manipulated to fit a long integer
    allnormalkeys = false; //lets the findkeys function know that there will be minikeys in the hashmaps from now on
    unsigned int miniseqlen = 10; //currently using 10 because why not.
    unsigned int iterations = (librarylength - (librarylength % miniseqlen)) / miniseqlen;
    unsigned int remainder = librarylength % miniseqlen;
    char miniseq[miniseqlen+1];
    char miniseqremainder[remainder+1];
    numerickey = 0;

    //takes the first 10 (or w/e miniseqlen is set to) nucleotides, converts to numeric (e.g. ATGAT = 12412), adds to numerickey (starts at 0). Repeats for next 10 ... until entire library finished.
    for (i = 0; i < iterations; i++){
      for (j = 0; j < miniseqlen; j++){
        if(library[i * miniseqlen + j] == 'A'){
          miniseq[j] = '1';
        } else if (library[i * miniseqlen + j] == 'T'){
          miniseq[j] = '2';
        } else if (library[i * miniseqlen + j] == 'C'){
          miniseq[j] = '3';
        } else if (library[i * miniseqlen + j] == 'G'){
          miniseq[j] = '4';
        } else { //currently sets all non ATCG to 5.
          miniseq[j] = '5';
          printf("Non-standard nucleotide detected in library.\n");
        }
      }
      miniseq[miniseqlen]='\0';
      numerickey = numerickey + strtoull(miniseq, NULL, 10);
    }
    for (i = 0; i < remainder; i++){
      if(library[iterations * miniseqlen + i] == 'A'){
        miniseqremainder[i] = '1';
      } else if (library[iterations * miniseqlen + i] == 'T'){
        miniseqremainder[i] = '2';
      } else if (library[iterations * miniseqlen + i] == 'C'){
        miniseqremainder[i] = '3';
      } else if (library[iterations * miniseqlen + i] == 'G'){
        miniseqremainder[i] = '4';
      } else { //currently sets all non ATCG to 5.
        miniseqremainder[i] = '5';
        printf("Non-standard nucleotide detected in library.\n");
      }
    }
    miniseqremainder[remainder]='\0';
    numerickey = numerickey + strtoull(miniseqremainder, NULL, 10);
  }
  return numerickey;
}

node *findnextnode(hashmap* restrict array, unsigned int index){ //finds the first non-null node starting from index on hashtable
  node *temp;
  for(unsigned int i=index; i<prime; i++){ //loops from array[index] to array[extern unsigned int prime]
    temp = array[i].head;
    if(temp){ //breaks loop if a non-null node was found, otherwise index++ and continues
      return temp;
    }
  }
  return NULL; //no more nodes
}

node *linkentirelist(hashmap* restrict array){ //given a hashtable of singly linked lists, link the whole thing and return the head node

  node *headnode, *temp, *current;

  headnode = findnextnode(array, 0);
  if(headnode == NULL){ //makes sure there was actually a headnode, otherwise exit this function
    printf("Error, no nodes were found in the specified hashtable.\n");
    return NULL;
  }
  unsigned int currentindex = headnode->key % prime; //finds the index of headnode and starts there
  temp = headnode;
  while(temp){ //This loop condition can be any condition that should never be reached. It is just set to (temp!=NULL) as a safety.
    while(temp){ //this loop sends current to last node of linked list for currentindex
      current = temp;
      temp = temp->next;
    }
    if (currentindex == prime){ //exits loop to prevent you from calling findnextnode on out-of-bounds index
      return headnode;
    }
    currentindex++;
    temp = findnextnode(array, currentindex); //set temp to next head node in hashtable
    if(!temp){ //prevents you from trying to call NULL->key
      return headnode;
    }
    currentindex = temp->key % prime;
    current->next = temp;
  }
    return headnode; //this should probably never be reached, but it's here just in case.
}

splitlist split_list(node* head){
  splitlist newlists;
  node *temp = head;
  node *prev = NULL;
  newlists.list1 = NULL;
  newlists.list2 = NULL;

  while(temp){
    prev = temp;
    temp = temp->next;
    if(prev->count == 1){ //all count = 1 goes into list 1
      if(newlists.list1){
        prev->next = newlists.list1;
        newlists.list1 = prev;
      } else {
        prev->next = NULL;
        newlists.list1 = prev;
      }
    } else { //all count != 1 goes into list 2
      if(newlists.list2){
        prev->next = newlists.list2;
        newlists.list2 = prev;
      } else {
        prev->next = NULL;
        newlists.list2 = prev;
      }
    }
  }
  return newlists;
}
