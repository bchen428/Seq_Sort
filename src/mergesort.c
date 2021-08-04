#include"main.h"
#include"structures.h"

//iterative merge sort
void merge(node** start1, node** end1, node** start2, node** end2){
  node *temp, *temp2;
  if ((*start1)->count < (*start2)->count){ //makes sure the first node of second list is higher
    temp2 = *start1;
    *start1 = *start2;
    *start2 = temp2;
    temp2 = *end1;
    *end1 = *end2;
    *end2 = temp2;
  }

  // Merging remaining nodes
  node *astart = *start1;
  node *aend = *end1;
  node *bstart = *start2;
  node *bendnext = (*end2)->next;
  while (astart != aend && bstart != bendnext) {
    if (astart->next->count < bstart->count){
      temp = bstart->next;
      bstart->next = astart->next;
      astart->next = bstart;
      bstart = temp;
    }
    astart = astart->next;
  }
  if (astart == aend){
    astart->next = bstart;
  } else {
    *end2 = *end1;
  }
}

void mergesort(node** head){ //sort linked list by changing pointers
  if (head == NULL){
    return;
  }
  node *start1 = NULL, *end1 = NULL, *start2 = NULL, *end2 = NULL, *prevend = NULL;

  node *current = *head;
  unsigned int len = 0;
  while (current != NULL) { //gets length of list
      current = current->next;
      len++;
  }

  for (unsigned int gap = 1; gap < len; gap += gap) {
    start1 = *head;
    while (start1){
      // If this is first iteration
      bool isFirstIter = 0;
      if (start1 == *head){
        isFirstIter = 1;
      }

      // First part for merging
      int counter = gap;
      end1 = start1;
      while (--counter && end1->next){
        end1 = end1->next;
      }

      // Second part for merging
      start2 = end1->next;
      if (!start2){
        break;
      }
      counter = gap;
      end2 = start2;
      while (--counter && end2->next){
        end2 = end2->next;
      }

      // To store for next iteration.
      node *temp = end2->next;

      // Merging two parts.
      merge(&start1, &end1, &start2, &end2);

      // Update head for first iteration, else
      // append after previous list
      if (isFirstIter){
        *head = start1;
      } else {
        prevend->next = start1;
      }
      prevend = end2;
      start1 = temp;
    }
    prevend->next = start1;
  }
}
