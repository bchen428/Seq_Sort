#ifndef STRUCTURES_HEADER
#define STRUCTURES_HEADER

typedef struct listNode {
  unsigned int count;
  unsigned int totalcount;
  struct listNode *next;
} listnode;

typedef struct Node {
  char *library;
  u_int64_t key;
  unsigned int count;
  struct Node *next;
} node;

typedef struct hashmap {
  node *head;
} hashmap;

typedef struct splitlist {
  node *list1;
  node *list2;
} splitlist;

#endif
