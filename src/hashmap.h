#ifndef HASHMAP_HEADER
#define HASHMAP_HEADER

int findkey(node* restrict list, u_int64_t key, char* restrict library);
node *get_node(node* restrict list, unsigned int index);
void insert(hashmap* restrict array, u_int64_t key, char* restrict library);
node *move(hashmap* restrict array, node* existingnode);
void initialize(hashmap* restrict array);
u_int64_t library_to_key(char* restrict library, unsigned int librarylength);
node *findnextnode(hashmap* restrict array, unsigned int index);
node *linkentirelist(hashmap* restrict array);
splitlist split_list(node* head);

#endif
