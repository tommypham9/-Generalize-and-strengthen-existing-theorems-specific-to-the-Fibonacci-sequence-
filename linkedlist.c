#include <stdio.h>
#include <stdlib.h>
#include "linkedlist.h"
#include <stdbool.h>
#include <gmp.h>
struct Node* createNode(mpz_t data){
    struct Node* newNode=(struct Node*)malloc(sizeof(struct Node));
    if(!newNode){
        printf("Memory allocation error\n");
        return NULL;
    }
    mpz_init(newNode->data);
    mpz_set(newNode->data,data);
    newNode->next=NULL;
    return newNode;
}

void insertAtBeginning(struct Node** head, mpz_t data){
    struct Node* newNode= createNode(data);
    newNode->next= *head;
    *head= newNode;

}
void displayList(struct Node* head){
    struct Node* current= head;
    if(current==NULL){
        printf("list is empty");
    }
    else{
    while(current!=NULL){
        gmp_printf("%Zd->", current->data);
        current=current->next;
    }
    printf("NULL\n");
    }
    
}

int listSize(struct Node* head){
    struct Node* current= head;
    int count=0;
    while(current!=NULL){
        count++;
        current=current->next;
    }
    return count;
}

int contains(struct Node* head, mpz_t data){
    struct Node* current = head;
    while (current != NULL){
        if(mpz_cmp(current->data,data)==0){
            return 1;
        }
        current= current->next;
    }
    return 0;
}

void deleteList(struct Node** head){
    struct Node* current=  *head;
    struct Node* next;
    while(current != NULL){
        next=current->next;
        mpz_clear(current->data);
        free(current);
        current=next;
    }
    *head=NULL;
}

void deleteNode(struct Node** head_ref, mpz_t key) {
    struct Node* temp = *head_ref;
    struct Node* prev = NULL;

    // If head node itself holds the key to be deleted
    if (temp != NULL && mpz_cmp(temp->data, key) == 0) {
        *head_ref = temp->next; // Change head
        mpz_clear(temp->data); // Free GMP data
        free(temp); // Free old head
        return;
    }

    // Search for the key to be deleted, keep track of the previous node
    while (temp != NULL && mpz_cmp(temp->data, key) != 0) {
        prev = temp;
        temp = temp->next;
    }

    // If key was not present in the linked list
    if (temp == NULL) return;

    // Unlink the node from the linked list
    prev->next = temp->next;

    // Free the memory allocated to the node
    mpz_clear(temp->data); // Free GMP data
    free(temp);
}