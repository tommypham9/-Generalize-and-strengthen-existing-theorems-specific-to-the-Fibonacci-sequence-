
#include <stdio.h>
#include <gmp.h> 
#include <stdlib.h>
#include "linkedlist.h"
#include <stdbool.h>
#include <math.h>
// Assume Node and related linked list functions (insertAtBeginning, contains, deleteList, etc.) are already defined.

void findDivisors(mpz_t num, struct Node** head) {
    mpz_t i;
    mpz_init(i);
    mpz_t sqrt_num;
    mpz_init(sqrt_num);
    mpz_sqrt(sqrt_num, num);

    for (mpz_set_ui(i, 1); mpz_cmp(i, sqrt_num) <= 0; mpz_add_ui(i, i, 1)) {
        if (mpz_divisible_p(num, i)) {
            insertAtBeginning(head, i);
            mpz_t div;
            mpz_init(div);
            mpz_divexact(div, num, i);
            if (mpz_cmp(i, div) != 0) {
                insertAtBeginning(head, div);
            }
            mpz_clear(div);
        }
    }

    mpz_clear(i);
    mpz_clear(sqrt_num);
}

/*void theNthLucasNum(int n, int p, int q, mpz_t result) {
    mpz_t L0, L1, L_temp;
    mpz_init_set_ui(L0, p);
    mpz_init_set_ui(L1, q);
    mpz_init(L_temp);

    if (n == 0) {
        mpz_set(result, L0);
    } else if (n == 1) {
        mpz_set(result, L1);
    } else {
        for (int i = 2; i <= n; ++i) {
            mpz_add(L_temp, L0, L1); //L_temp= L0 + L1
            mpz_set(L0, L1); //L0= L1
            mpz_set(L1, L_temp); // L1= L_temp
        }
        mpz_set(result, L1);
    }

    mpz_clear(L0);
    mpz_clear(L1);
    mpz_clear(L_temp);
}*/


// Function to multiply two 2x2 matrices
void matrix_multiply(mpz_t result[2][2], mpz_t mat1[2][2], mpz_t mat2[2][2]) {
    mpz_t a, b, c, d;
    mpz_init(a);
    mpz_init(b);
    mpz_init(c);
    mpz_init(d);

    mpz_mul(a, mat1[0][0], mat2[0][0]);
    mpz_addmul(a, mat1[0][1], mat2[1][0]);
    
    mpz_mul(b, mat1[0][0], mat2[0][1]);
    mpz_addmul(b, mat1[0][1], mat2[1][1]);

    mpz_mul(c, mat1[1][0], mat2[0][0]);
    mpz_addmul(c, mat1[1][1], mat2[1][0]);
    
    mpz_mul(d, mat1[1][0], mat2[0][1]);
    mpz_addmul(d, mat1[1][1], mat2[1][1]);

    mpz_set(result[0][0], a);
    mpz_set(result[0][1], b);
    mpz_set(result[1][0], c);
    mpz_set(result[1][1], d);

    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(c);
    mpz_clear(d);
}

// Function to exponentiate a 2x2 matrix to the power n
void matrix_power(mpz_t result[2][2], mpz_t base[2][2], unsigned long long n) {
    mpz_t temp[2][2];
    mpz_init(temp[0][0]);
    mpz_init(temp[0][1]);
    mpz_init(temp[1][0]);
    mpz_init(temp[1][1]);

    // Initialize result as the identity matrix
    mpz_set_ui(result[0][0], 1);
    mpz_set_ui(result[0][1], 0);
    mpz_set_ui(result[1][0], 0);
    mpz_set_ui(result[1][1], 1);

    while (n > 0) {
        if (n % 2 == 1) {
            matrix_multiply(temp, result, base);
            mpz_set(result[0][0], temp[0][0]);
            mpz_set(result[0][1], temp[0][1]);
            mpz_set(result[1][0], temp[1][0]);
            mpz_set(result[1][1], temp[1][1]);
        }
        matrix_multiply(temp, base, base);
        mpz_set(base[0][0], temp[0][0]);
        mpz_set(base[0][1], temp[0][1]);
        mpz_set(base[1][0], temp[1][0]);
        mpz_set(base[1][1], temp[1][1]);
        n /= 2;
    }

    mpz_clear(temp[0][0]);
    mpz_clear(temp[0][1]);
    mpz_clear(temp[1][0]);
    mpz_clear(temp[1][1]);
}

void generalizedNthFibNum(unsigned long long n, mpz_t a, mpz_t b, mpz_t result) {
    if (n == 0) {
        mpz_set(result, a);
        return;
    }
    if (n == 1) {
        mpz_set(result, b);
        return;
    }

    mpz_t base[2][2];
    mpz_init_set_ui(base[0][0], 1);
    mpz_init_set_ui(base[0][1], 1);
    mpz_init_set_ui(base[1][0], 1);
    mpz_init_set_ui(base[1][1], 0);

    mpz_t result_matrix[2][2];
    mpz_init(result_matrix[0][0]);
    mpz_init(result_matrix[0][1]);
    mpz_init(result_matrix[1][0]);
    mpz_init(result_matrix[1][1]);

    matrix_power(result_matrix, base, n-1);

    mpz_t temp_a, temp_b;
    mpz_init(temp_a);
    mpz_init(temp_b);

    // result = result_matrix[0][0] * b + result_matrix[0][1] * a
    mpz_mul(temp_a, result_matrix[0][0], b);
    mpz_mul(temp_b, result_matrix[0][1], a);
    mpz_add(result, temp_a, temp_b);

    mpz_clear(temp_a);
    mpz_clear(temp_b);
    mpz_clear(base[0][0]);
    mpz_clear(base[0][1]);
    mpz_clear(base[1][0]);
    mpz_clear(base[1][1]);
    mpz_clear(result_matrix[0][0]);
    mpz_clear(result_matrix[0][1]);
    mpz_clear(result_matrix[1][0]);
    mpz_clear(result_matrix[1][1]);
}

void gcd(const mpz_t a, const mpz_t b, mpz_t result) {
    mpz_t temp_a, temp_b, temp;
    mpz_init_set(temp_a, a); // temp_a = a
    mpz_init_set(temp_b, b); // temp_b = b
    mpz_init(temp);

    while (mpz_cmp_ui(temp_b, 0) != 0) {
        mpz_set(temp, temp_b);
        mpz_mod(temp_b, temp_a, temp_b);
        mpz_set(temp_a, temp);
    }
    
    if (mpz_cmp_ui(temp_a, 0) < 0) {
        mpz_abs(temp_a, temp_a);
    }
    
    mpz_set(result, temp_a);

    mpz_clear(temp);
    mpz_clear(temp_a);
    mpz_clear(temp_b);
}


int isPrime(const mpz_t n){
    mpz_t sqrt_num;
    mpz_init(sqrt_num);
    mpz_sqrt(sqrt_num,n);
    mpz_t i;
    mpz_init(i);
    if(mpz_cmp_ui(n,1)<=0){
        return 0;
    }
    if(mpz_cmp_ui(n,3)<=0){
        return 1;
    }
    for(mpz_set_ui(i,2);mpz_cmp(i,sqrt_num)<=0; mpz_add_ui(i,i,1)){
        if(mpz_divisible_p(n,i)!=0){
            return 0;
        }
    }
    mpz_clear(sqrt_num);
    mpz_clear(i);
    return 1;
}

struct Node* primeFactors(const mpz_t original_n) {
    struct Node* div_list = NULL;
    mpz_t n, temp, i, sqrt_num, temp1;
    mpz_inits(n, temp, i, sqrt_num, temp1, NULL);

    mpz_set(n, original_n); // Use a copy of original_n
    mpz_set_ui(temp, 2);
    mpz_sqrt(sqrt_num, n);

    while (mpz_divisible_ui_p(n, 2) != 0) {
        if (contains(div_list, temp) == 0) {
            insertAtBeginning(&div_list, temp);
        }
        mpz_div_ui(n, n, 2);
    }

    mpz_set_ui(i, 3);
    mpz_mul(temp1, i, i);
    while (mpz_cmp(temp1, n) <= 0) {
        while (mpz_divisible_p(n, i) != 0) {
            if (contains(div_list, i) == 0) {
                insertAtBeginning(&div_list, i);
            }
            mpz_divexact(n, n, i);
        }
        mpz_add_ui(i, i, 2);
        mpz_mul(temp1, i, i);
    }

    if (mpz_cmp_ui(n, 2) > 0 && contains(div_list, n) == 0) {
        insertAtBeginning(&div_list, n);
    }

    mpz_clears(n, temp, i, sqrt_num, temp1, NULL);
    return div_list;
}

struct Node* k_1Moduli(unsigned long long h_period){
    unsigned long long k_1= h_period;
    struct Node* m1_head=NULL;
    struct Node* f1_head=NULL;
    struct Node* div1_head=NULL;
    struct Node* head=NULL;
    mpz_t a,b,c;                //remember to free those 
    mpz_init_set_ui(a,2);
    mpz_init_set_ui(b,1);
    mpz_init_set_ui(c,0);
    if(k_1%4==2 && k_1 >=6){
        mpz_t lucas_number;
        mpz_init(lucas_number);
        generalizedNthFibNum(k_1/2,a,b,lucas_number );
        findDivisors(lucas_number,&m1_head);
        insertAtBeginning(&div1_head,b);
        for(unsigned long long i=2;i*i<=k_1;i++){
            if(k_1%i==0){
                mpz_t convert_i;
                mpz_init_set_ui(convert_i,i);
                insertAtBeginning(&div1_head,convert_i);
                if(i!=k_1/i){
                    mpz_t convert_j;
                    mpz_init_set_ui(convert_j,k_1/i);
                    insertAtBeginning(&div1_head, convert_j);
                    mpz_clear(convert_j);
                }
                mpz_clear(convert_i);
        }
    }
   
        struct Node* current1= div1_head;
        while(current1!=NULL){
            mpz_t result;
            mpz_init(result);
            unsigned long long data= mpz_get_ui(current1->data);
            generalizedNthFibNum(data,c,b,result);
            insertAtBeginning(&f1_head,result);         //f list
            current1=current1->next;
            mpz_clear(result);
            }
        struct Node* m1_current=m1_head;
        struct Node* f1_current= f1_head;
        while(m1_current!=NULL){        //problem right here
        int count=0;
        f1_current=f1_head;
        while(f1_current!=NULL && count==0){
            if(mpz_divisible_p(f1_current->data, m1_current->data)==0){
                f1_current=f1_current->next;
            }
            else{
                count++;
            }
        }
            if(count==0 && contains(head, m1_current->data)==0){
                insertAtBeginning(&head, m1_current->data);
            }
            m1_current= m1_current->next;
        }
        mpz_clear(lucas_number);
    }
    else if(k_1%8 ==4 && k_1>=6){
        mpz_t fib_num, lucas_num;       //remember to free those later- i did ! 
        mpz_init(fib_num);
        mpz_init(lucas_num);
        generalizedNthFibNum(k_1/2,c,b,fib_num);
        generalizedNthFibNum(k_1/4,a,b,lucas_num);
        findDivisors(fib_num,&m1_head);
        
        struct Node* check_m= m1_head;
        struct Node* next_node;
        while(check_m!=NULL){
            next_node=check_m->next;
            if(mpz_divisible_p(lucas_num,check_m->data)!=0){
                deleteNode(&m1_head,check_m->data);
            }
            check_m=next_node;
        }
        insertAtBeginning(&div1_head,b);
        for(unsigned long long i=2;i*i<=k_1/2;i++){
            if((k_1/2)%i==0 && k_1/4 != i){
                mpz_t convert_i;
                mpz_init_set_ui(convert_i,i);
                insertAtBeginning(&div1_head,convert_i);
                if(i!=(k_1/2)/i && (k_1/2)/i != k_1/4){
                    mpz_t convert_j;
                    mpz_init_set_ui(convert_j,k_1/i);
                    insertAtBeginning(&div1_head, convert_j);
                    mpz_clear(convert_j);
                }
                mpz_clear(convert_i);
        }
        }
        struct Node* current_div= div1_head;
        while(current_div!=NULL){
            mpz_t result;
            mpz_init(result);
            unsigned long long data= mpz_get_ui(current_div->data);
            generalizedNthFibNum(data,c,b,result);
            insertAtBeginning(&f1_head,result); //f list
            current_div=current_div->next;
            mpz_clear(result);
        }
        struct Node* m1_current= m1_head;
        while(m1_current!=NULL){
            int count=0;
            struct Node* f1_current= f1_head;
            while(f1_current!=NULL && count ==0){
                if(mpz_divisible_p(f1_current->data,m1_current->data)==0){
                    f1_current=f1_current->next;
                }
                else{
                    count++;
                }
            }
            if(count==0 && contains(head,m1_current->data)==0){
                insertAtBeginning(&head,m1_current->data);
            }
            m1_current=m1_current->next;
        }
        mpz_clear(fib_num);
        mpz_clear(lucas_num);
        }
        else{
            mpz_t fib_num;
            mpz_init(fib_num);
            generalizedNthFibNum(k_1/2,c,b,fib_num);
            findDivisors(fib_num,&m1_head);     //m list
            insertAtBeginning(&div1_head, b);
            for(unsigned long long i=2; i*i<=k_1/2;i++){
                if((k_1/2)%i==0){
                mpz_t convert_i;
                mpz_init_set_ui(convert_i,i);
                insertAtBeginning(&div1_head,convert_i);
                    if(i!=(k_1/2)/i){
                        mpz_t convert_j;
                        mpz_init_set_ui(convert_j,(k_1/2)/i);
                        insertAtBeginning(&div1_head, convert_j);
                        mpz_clear(convert_j);
                }
                mpz_clear(convert_i);
        }
            }
            struct Node* current_div= div1_head;
            while(current_div!=NULL){
                mpz_t result;
                mpz_init(result);
                unsigned long long data= mpz_get_ui(current_div->data);
                generalizedNthFibNum(data,c,b,result);
                insertAtBeginning(&f1_head,result);         // f list
                current_div=current_div->next;
                mpz_clear(result);
            }
            struct Node* m1_current= m1_head;
            while(m1_current!=NULL){
                int count=0;
                struct Node* f1_current=f1_head;
                while(f1_current!=NULL && count==0){
                    if(mpz_divisible_p(f1_current->data, m1_current->data)==0){
                        f1_current=f1_current->next;
                    }
                    else{
                        count++;
                    }
                }
                
                if(count==0 && contains(head, m1_current->data)==0){
                    insertAtBeginning(&head,m1_current->data);
                }
                m1_current=m1_current->next;
            }
            mpz_clear(fib_num);
        }
    deleteList(&m1_head);
    deleteList(&f1_head);
    deleteList(&div1_head);
    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(c);
    return head;
}

struct Node* potentialModuli(unsigned long long h_period,unsigned long long det){
    unsigned long long k_2, k_3, k_4;
    k_2=2*h_period;
    k_3=5*h_period;
    k_4=10*h_period;
    struct Node* head=NULL;
    mpz_t a,b,c;                //remember to free those 
    mpz_init_set_ui(a,2);
    mpz_init_set_ui(b,1);
    mpz_init_set_ui(c,0);
    if(h_period%2==1){
        struct Node* m1_head= NULL;
        struct Node* f1_head= NULL;
        struct Node* div1_head=NULL;
        if(k_2 %4 ==2 && k_2>=6){
            mpz_t lucas_number;
            mpz_init(lucas_number);
            generalizedNthFibNum(k_2/2,a,b,lucas_number );
            findDivisors(lucas_number,&m1_head);
            insertAtBeginning(&div1_head,b);
            for(unsigned long long i=2;i*i<=k_2;i++){
                if(k_2%i==0){
                    mpz_t convert_i;
                    mpz_init_set_ui(convert_i,i);
                    insertAtBeginning(&div1_head,convert_i);
                    if(i!=k_2/i){
                        mpz_t convert_j;
                        mpz_init_set_ui(convert_j,k_2/i);
                        insertAtBeginning(&div1_head, convert_j);
                        mpz_clear(convert_j);
                }
                mpz_clear(convert_i);
        }
    }
   
        struct Node* current1= div1_head;
        while(current1!=NULL){
            mpz_t result;
            mpz_init(result);
            unsigned long long data= mpz_get_ui(current1->data);
            generalizedNthFibNum(data,c,b,result);
            insertAtBeginning(&f1_head,result);         //f list
            current1=current1->next;
            mpz_clear(result);
            }
        struct Node* m1_current=m1_head;
        struct Node* f1_current= f1_head;
        while(m1_current!=NULL){        //problem right here
        int count=0;
        f1_current=f1_head;
        while(f1_current!=NULL && count==0){
            if(mpz_divisible_p(f1_current->data, m1_current->data)==0){
                f1_current=f1_current->next;
            }
            else{
                count++;
            }
        }
            if(count==0 && contains(head, m1_current->data)==0){
                insertAtBeginning(&head, m1_current->data);
            }
            m1_current= m1_current->next;
        }
        mpz_clear(lucas_number);
        }
        deleteList(&m1_head);
        deleteList(&f1_head);
        deleteList(&div1_head);

        //break right here

        if(k_4%4==2 && k_4>=6){
            mpz_t lucas_number;
            mpz_init(lucas_number);
            generalizedNthFibNum(k_4/2,a,b,lucas_number );
            findDivisors(lucas_number,&m1_head);
            insertAtBeginning(&div1_head,b);
            for(unsigned long long i=2;i*i<=k_4;i++){
                if(k_4%i==0){   //k_2
                    mpz_t convert_i;
                    mpz_init_set_ui(convert_i,i);
                    insertAtBeginning(&div1_head,convert_i);
                    if(i!=k_4/i){
                        mpz_t convert_j;
                        mpz_init_set_ui(convert_j,k_4/i);
                        insertAtBeginning(&div1_head, convert_j);
                        mpz_clear(convert_j);
                }
                mpz_clear(convert_i);
        }
    }
   
        struct Node* current1= div1_head;
        while(current1!=NULL){
            mpz_t result;
            mpz_init(result);
            unsigned long long data= mpz_get_ui(current1->data);
            generalizedNthFibNum(data,c,b,result);
            insertAtBeginning(&f1_head,result);         //f list
            current1=current1->next;
            mpz_clear(result);
            }
        struct Node* m1_current=m1_head;
        struct Node* f1_current= f1_head;
        while(m1_current!=NULL){        //problem right here
        int count=0;
        f1_current=f1_head;
        while(f1_current!=NULL && count==0){
            if(mpz_divisible_p(f1_current->data, m1_current->data)==0){
                f1_current=f1_current->next;
            }
            else{
                count++;
            }
        }
            if(count==0 && contains(head, m1_current->data)==0){
                insertAtBeginning(&head, m1_current->data);
            }
            m1_current= m1_current->next;
        }
        mpz_clear(lucas_number);
        }
        deleteList(&m1_head);
        deleteList(&f1_head);       //free memory 
        deleteList(&div1_head);     

    }
    else{
        struct Node* m1_head=NULL;
        struct Node* f1_head=NULL;
        struct Node* div1_head=NULL;

        if(k_2%4==2 && k_2 >=6){
            mpz_t lucas_number;
            mpz_init(lucas_number);
            generalizedNthFibNum(k_2/2,a,b,lucas_number );
            findDivisors(lucas_number,&m1_head);
            insertAtBeginning(&div1_head,b);
            for(unsigned long long i=2;i*i<=k_2;i++){
                if(k_2%i==0){
                    mpz_t convert_i;
                    mpz_init_set_ui(convert_i,i);
                    insertAtBeginning(&div1_head,convert_i);
                    if(i!=k_2/i){
                        mpz_t convert_j;
                        mpz_init_set_ui(convert_j,k_2/i);
                        insertAtBeginning(&div1_head, convert_j);
                        mpz_clear(convert_j);
                }
                mpz_clear(convert_i);
        }
    }
   
        struct Node* current1= div1_head;
        while(current1!=NULL){
            mpz_t result;
            mpz_init(result);
            unsigned long long data= mpz_get_ui(current1->data);
            generalizedNthFibNum(data,c,b,result);
            insertAtBeginning(&f1_head,result);         //f list
            current1=current1->next;
            mpz_clear(result);
            }
        
        struct Node* m1_current=m1_head;
        struct Node* f1_current= f1_head;
        while(m1_current!=NULL){        //problem right here
        int count=0;
        f1_current=f1_head;
        while(f1_current!=NULL && count==0){
            if(mpz_divisible_p(f1_current->data, m1_current->data)==0){
                f1_current=f1_current->next;
            }
            else{
                count++;
            }
        }
            if(count==0 && contains(head, m1_current->data)==0){
                insertAtBeginning(&head, m1_current->data);
            }
            m1_current= m1_current->next;
        }
        mpz_clear(lucas_number);
        }

        else if(k_2 %8==4 && k_2>=6){
            mpz_t fib_num, lucas_num;       //remember to free those later- i did ! 
            mpz_init(fib_num);
            mpz_init(lucas_num);
            generalizedNthFibNum(k_2/2,c,b,fib_num);
            generalizedNthFibNum(k_2/4,a,b,lucas_num);
            findDivisors(fib_num,&m1_head);
        //k_1
            struct Node* check_m= m1_head;
            struct Node* next_node;
            while(check_m!=NULL){
            next_node=check_m->next;
            if(mpz_divisible_p(lucas_num,check_m->data)!=0){
                deleteNode(&m1_head,check_m->data);
            }
            check_m=next_node;
        }
        insertAtBeginning(&div1_head,b);
        for(unsigned long long i=2;i*i<=k_2/2;i++){
            if((k_2/2)%i==0 && k_2/4 != i){
                mpz_t convert_i;
                mpz_init_set_ui(convert_i,i);
                insertAtBeginning(&div1_head,convert_i);
                if(i!=(k_2/2)/i && (k_2/2)/i != k_2/4){
                    mpz_t convert_j;
                    mpz_init_set_ui(convert_j,k_2/i);
                    insertAtBeginning(&div1_head, convert_j);
                    mpz_clear(convert_j);
                }
                mpz_clear(convert_i);
        }
        }
        struct Node* current_div= div1_head;
        while(current_div!=NULL){
            mpz_t result;
            mpz_init(result);
            unsigned long long data= mpz_get_ui(current_div->data);
            generalizedNthFibNum(data,c,b,result);
            insertAtBeginning(&f1_head,result); //f list
            current_div=current_div->next;
            mpz_clear(result);
        }
        struct Node* m1_current= m1_head;
        while(m1_current!=NULL){
            int count=0;
            struct Node* f1_current= f1_head;
            while(f1_current!=NULL && count ==0){
                if(mpz_divisible_p(f1_current->data,m1_current->data)==0){
                    f1_current=f1_current->next;
                }
                else{
                    count++;
                }
            }
            if(count==0 && contains(head,m1_current->data)==0){
                insertAtBeginning(&head,m1_current->data);
            }
            m1_current=m1_current->next;
        }
        mpz_clear(fib_num);
        mpz_clear(lucas_num);
        }
        else{
            mpz_t fib_num;
            mpz_init(fib_num);
            generalizedNthFibNum(k_2/2,c,b,fib_num);        //k_1
            findDivisors(fib_num,&m1_head);     //m list
            insertAtBeginning(&div1_head, b);
            for(unsigned long long i=2; i*i<=k_2/2;i++){
                if((k_2/2)%i==0){
                mpz_t convert_i;
                mpz_init_set_ui(convert_i,i);
                insertAtBeginning(&div1_head,convert_i);
                    if(i!=(k_2/2)/i){
                        mpz_t convert_j;
                        mpz_init_set_ui(convert_j,(k_2/2)/i);
                        insertAtBeginning(&div1_head, convert_j);
                        mpz_clear(convert_j);
                }
                mpz_clear(convert_i);
        }
            }
            struct Node* current_div= div1_head;
            while(current_div!=NULL){
                mpz_t result;
                mpz_init(result);
                unsigned long long data= mpz_get_ui(current_div->data);
                generalizedNthFibNum(data,c,b,result);
                insertAtBeginning(&f1_head,result);         // f list
                current_div=current_div->next;
                mpz_clear(result);
            }
            struct Node* m1_current= m1_head;
            while(m1_current!=NULL){
                int count=0;
                struct Node* f1_current=f1_head;
                while(f1_current!=NULL && count==0){
                    if(mpz_divisible_p(f1_current->data, m1_current->data)==0){
                        f1_current=f1_current->next;
                    }
                    else{
                        count++;
                    }
                }
                
                if(count==0 && contains(head, m1_current->data)==0){
                    insertAtBeginning(&head,m1_current->data);
                }
                m1_current=m1_current->next;
            }
            mpz_clear(fib_num);
        }
    deleteList(&m1_head);
    deleteList(&f1_head);       //free memory 
    deleteList(&div1_head);

    if(k_3%4==2 && k_3>=6){
        mpz_t lucas_number;     
            mpz_init(lucas_number);
            generalizedNthFibNum(k_3/2,a,b,lucas_number );
            findDivisors(lucas_number,&m1_head);
            insertAtBeginning(&div1_head,b);
            for(unsigned long long i=2;i*i<=k_3;i++){
                if(k_3%i==0){
                    mpz_t convert_i;
                    mpz_init_set_ui(convert_i,i);
                    insertAtBeginning(&div1_head,convert_i);
                    if(i!=k_3/i){
                        mpz_t convert_j;
                        mpz_init_set_ui(convert_j,k_3/i);
                        insertAtBeginning(&div1_head, convert_j);       //k_2
                        mpz_clear(convert_j);
                }
                mpz_clear(convert_i);
        }
    }
   
        struct Node* current1= div1_head;
        while(current1!=NULL){
            mpz_t result;
            mpz_init(result);
            unsigned long long data= mpz_get_ui(current1->data);
            generalizedNthFibNum(data,c,b,result);
            insertAtBeginning(&f1_head,result);         //f list
            current1=current1->next;
            mpz_clear(result);
            }
        struct Node* m1_current=m1_head;
        struct Node* f1_current= f1_head;
        while(m1_current!=NULL){        //problem right here
        int count=0;
        f1_current=f1_head;
        while(f1_current!=NULL && count==0){
            if(mpz_divisible_p(f1_current->data, m1_current->data)==0){
                f1_current=f1_current->next;
            }
            else{
                count++;
            }
        }
            if(count==0 && contains(head, m1_current->data)==0){
                insertAtBeginning(&head, m1_current->data);
            }
            m1_current= m1_current->next;
        }
        mpz_clear(lucas_number);
    }
    else if(k_3%8 ==4 && k_3>=6){
            mpz_t fib_num, lucas_num;       //remember to free those later- i did ! 
            mpz_init(fib_num);
            mpz_init(lucas_num);
            generalizedNthFibNum(k_3/2,c,b,fib_num);
            generalizedNthFibNum(k_3/4,a,b,lucas_num);
            findDivisors(fib_num,&m1_head);
        //k_1
            struct Node* check_m= m1_head;
            struct Node* next_node;
            while(check_m!=NULL){
            next_node=check_m->next;
            if(mpz_divisible_p(lucas_num,check_m->data)!=0){
                deleteNode(&m1_head,check_m->data);
            }
            check_m=next_node;
        }
        insertAtBeginning(&div1_head,b);
        for(unsigned long long i=2;i*i<=k_3/2;i++){
            if((k_3/2)%i==0 && k_3/4 != i){
                mpz_t convert_i;
                mpz_init_set_ui(convert_i,i);
                insertAtBeginning(&div1_head,convert_i);
                if(i!=(k_3/2)/i && (k_3/2)/i != k_3/4){     //k_2
                    mpz_t convert_j;
                    mpz_init_set_ui(convert_j,k_3/i);
                    insertAtBeginning(&div1_head, convert_j);
                    mpz_clear(convert_j);
                }
                mpz_clear(convert_i);
        }
        }
        struct Node* current_div= div1_head;
        while(current_div!=NULL){
            mpz_t result;
            mpz_init(result);
            unsigned long long data= mpz_get_ui(current_div->data);
            generalizedNthFibNum(data,c,b,result);
            insertAtBeginning(&f1_head,result); //f list
            current_div=current_div->next;
            mpz_clear(result);
        }
        struct Node* m1_current= m1_head;
        while(m1_current!=NULL){
            int count=0;
            struct Node* f1_current= f1_head;
            while(f1_current!=NULL && count ==0){
                if(mpz_divisible_p(f1_current->data,m1_current->data)==0){
                    f1_current=f1_current->next;
                }
                else{
                    count++;
                }
            }
            if(count==0 && contains(head,m1_current->data)==0){
                insertAtBeginning(&head,m1_current->data);
            }
            m1_current=m1_current->next;
        }
        mpz_clear(fib_num);
        mpz_clear(lucas_num);
    }
    else{
            mpz_t fib_num;
            mpz_init(fib_num);
            generalizedNthFibNum(k_3/2,c,b,fib_num);        //k_2
            findDivisors(fib_num,&m1_head);     //m list
            insertAtBeginning(&div1_head, b);
            for(unsigned long long i=2; i*i<=k_3/2;i++){
                if((k_3/2)%i==0){
                mpz_t convert_i;
                mpz_init_set_ui(convert_i,i);
                insertAtBeginning(&div1_head,convert_i);
                    if(i!=(k_3/2)/i){
                        mpz_t convert_j;
                        mpz_init_set_ui(convert_j,(k_3/2)/i);
                        insertAtBeginning(&div1_head, convert_j);
                        mpz_clear(convert_j);
                }
                mpz_clear(convert_i);
        }
            }
            struct Node* current_div= div1_head;
            while(current_div!=NULL){
                mpz_t result;
                mpz_init(result);
                unsigned long long data= mpz_get_ui(current_div->data);
                generalizedNthFibNum(data,c,b,result);
                insertAtBeginning(&f1_head,result);         // f list
                current_div=current_div->next;
                mpz_clear(result);
            }           //k_2
            struct Node* m1_current= m1_head;
            while(m1_current!=NULL){
                int count=0;
                struct Node* f1_current=f1_head;
                while(f1_current!=NULL && count==0){
                    if(mpz_divisible_p(f1_current->data, m1_current->data)==0){
                        f1_current=f1_current->next;
                    }
                    else{
                        count++;
                    }
                }
                
                if(count==0 && contains(head, m1_current->data)==0){
                    insertAtBeginning(&head,m1_current->data);
                }
                m1_current=m1_current->next;
            }
            mpz_clear(fib_num);
    }

    deleteList(&m1_head);
    deleteList(&f1_head);       //free memory 
    deleteList(&div1_head);

    if(k_4%4==2 && k_4>=6){
        mpz_t lucas_number;     
            mpz_init(lucas_number);
            generalizedNthFibNum(k_4/2,a,b,lucas_number );
            findDivisors(lucas_number,&m1_head);
            insertAtBeginning(&div1_head,b);
            for(unsigned long long i=2;i*i<=k_4;i++){
                if(k_4%i==0){
                    mpz_t convert_i;
                    mpz_init_set_ui(convert_i,i);
                    insertAtBeginning(&div1_head,convert_i);
                    if(i!=k_4/i){
                        mpz_t convert_j;
                        mpz_init_set_ui(convert_j,k_4/i);
                        insertAtBeginning(&div1_head, convert_j);       //k_3
                        mpz_clear(convert_j);
                }
                mpz_clear(convert_i);
        }
    }
   
        struct Node* current1= div1_head;
        while(current1!=NULL){
            mpz_t result;
            mpz_init(result);
            unsigned long long data= mpz_get_ui(current1->data);
            generalizedNthFibNum(data,c,b,result);
            insertAtBeginning(&f1_head,result);         //f list
            current1=current1->next;
            mpz_clear(result);
            }
        struct Node* m1_current=m1_head;
        struct Node* f1_current= f1_head;
        while(m1_current!=NULL){        //problem right here
        int count=0;
        f1_current=f1_head;
        while(f1_current!=NULL && count==0){
            if(mpz_divisible_p(f1_current->data, m1_current->data)==0){
                f1_current=f1_current->next;
            }
            else{
                count++;
            }
        }
            if(count==0 && contains(head, m1_current->data)==0){
                insertAtBeginning(&head, m1_current->data);
            }
            m1_current= m1_current->next;
        }
        mpz_clear(lucas_number);
    }
    else if(k_4%8==4 && k_4>=6){
        mpz_t fib_num, lucas_num;       //remember to free those later- i did ! 
            mpz_init(fib_num);
            mpz_init(lucas_num);
            generalizedNthFibNum(k_4/2,c,b,fib_num);
            generalizedNthFibNum(k_4/4,a,b,lucas_num);
            findDivisors(fib_num,&m1_head);
        //k_1
            struct Node* check_m= m1_head;
            struct Node* next_node;
            while(check_m!=NULL){
            next_node=check_m->next;
            if(mpz_divisible_p(lucas_num,check_m->data)!=0){
                deleteNode(&m1_head,check_m->data);
            }
            check_m=next_node;
        }
        insertAtBeginning(&div1_head,b);
        for(unsigned long long i=2;i*i<=k_4/2;i++){
            if((k_4/2)%i==0 && k_4/4 != i){
                mpz_t convert_i;
                mpz_init_set_ui(convert_i,i);
                insertAtBeginning(&div1_head,convert_i);
                if(i!=(k_4/2)/i && (k_4/2)/i != k_4/4){     //k_3
                    mpz_t convert_j;
                    mpz_init_set_ui(convert_j,k_4/i);
                    insertAtBeginning(&div1_head, convert_j);
                    mpz_clear(convert_j);
                }
                mpz_clear(convert_i);
        }
        }
        struct Node* current_div= div1_head;
        while(current_div!=NULL){
            mpz_t result;
            mpz_init(result);
            unsigned long long data= mpz_get_ui(current_div->data);
            generalizedNthFibNum(data,c,b,result);
            insertAtBeginning(&f1_head,result); //f list
            current_div=current_div->next;
            mpz_clear(result);
        }
        struct Node* m1_current= m1_head;
        while(m1_current!=NULL){
            int count=0;
            struct Node* f1_current= f1_head;
            while(f1_current!=NULL && count ==0){
                if(mpz_divisible_p(f1_current->data,m1_current->data)==0){
                    f1_current=f1_current->next;
                }
                else{
                    count++;
                }
            }
            if(count==0 && contains(head,m1_current->data)==0){
                insertAtBeginning(&head,m1_current->data);
            }
            m1_current=m1_current->next;
        }
        mpz_clear(fib_num);
        mpz_clear(lucas_num);
    }
    else{
            mpz_t fib_num;
            mpz_init(fib_num);
            generalizedNthFibNum(k_4/2,c,b,fib_num);        //k_3
            findDivisors(fib_num,&m1_head);     //m list
            insertAtBeginning(&div1_head, b);
            for(unsigned long long i=2; i*i<=k_4/2;i++){
                if((k_4/2)%i==0){
                mpz_t convert_i;
                mpz_init_set_ui(convert_i,i);
                insertAtBeginning(&div1_head,convert_i);
                    if(i!=(k_4/2)/i){
                        mpz_t convert_j;
                        mpz_init_set_ui(convert_j,(k_4/2)/i);
                        insertAtBeginning(&div1_head, convert_j);
                        mpz_clear(convert_j);
                }
                mpz_clear(convert_i);
        }
            }
            struct Node* current_div= div1_head;
            while(current_div!=NULL){
                mpz_t result;
                mpz_init(result);
                unsigned long long data= mpz_get_ui(current_div->data);
                generalizedNthFibNum(data,c,b,result);
                insertAtBeginning(&f1_head,result);         // f list
                current_div=current_div->next;
                mpz_clear(result);
            }           //k_2
            struct Node* m1_current= m1_head;
            while(m1_current!=NULL){
                int count=0;
                struct Node* f1_current=f1_head;
                while(f1_current!=NULL && count==0){
                    if(mpz_divisible_p(f1_current->data, m1_current->data)==0){
                        f1_current=f1_current->next;
                    }
                    else{
                        count++;
                    }
                }
                
                if(count==0 && contains(head, m1_current->data)==0){
                    insertAtBeginning(&head,m1_current->data);
                }
                m1_current=m1_current->next;
            }
            mpz_clear(fib_num);
    }
    deleteList(&m1_head);
    deleteList(&f1_head);       //free memory 
    deleteList(&div1_head);

    }
    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(c);
    return head;

}
struct Node* CalculateGenFibMod(struct Node* k1_list, struct Node* otherModuliList, mpz_t G_0,mpz_t G_1,unsigned long long period, mpz_t det){
    mpz_t fib_num, second_fib_num, converted_period;
    mpz_init(fib_num);
    mpz_init(second_fib_num);
    mpz_init(converted_period);
    mpz_set_ui(converted_period,period);
   // gmp_printf("\nconverted period is %Zd!!!!!!\n", converted_period);
    generalizedNthFibNum(period,G_0,G_1,fib_num);
    generalizedNthFibNum(period+1,G_0,G_1,second_fib_num);
    struct Node* final_list=NULL;
    struct Node* divisor= primeFactors(converted_period);
   // gmp_printf("\nafter converted period is %Zd!!!!!!\n", converted_period);
    struct Node* current_div= divisor;
    struct Node* current_k1=k1_list;
    struct Node* current_otherModuli= otherModuliList;
    while(current_k1!=NULL){
        int count=0;
        current_div=divisor;
        mpz_t result_gcd;       //free this later- i did 
        mpz_init(result_gcd);
        gcd(det, current_k1->data,result_gcd);
        if(mpz_cmp_ui(result_gcd,1)==0){
            insertAtBeginning(&final_list,current_k1->data);
        }
        else{
            while(current_div!=NULL && count==0){
                mpz_t G_d,G_d_1, temp_res, mod_result1, mod_result2;
                mpz_init(G_d);
                mpz_init(G_d_1);
                mpz_init(temp_res);
                mpz_init(mod_result1);
                mpz_init(mod_result2);
                mpz_div(temp_res, converted_period,current_div->data);
               // gmp_printf("\ntemp_res is %Zd\n converted period is %Zd\n", temp_res, converted_period);
                unsigned long long data= mpz_get_ui(temp_res);
                generalizedNthFibNum(data,G_0,G_1,G_d);
                generalizedNthFibNum(data+1,G_0,G_1,G_d_1);
               // gmp_printf("\nG_d is %Zd \nG_d_1 is %Zd\n G_0 is %Zd\n G_1 is %Zd\n data is %llu", G_d,G_d_1,G_0,G_1,data);
                mpz_mod(mod_result1,G_d,current_k1->data);
                mpz_mod(mod_result2,G_d_1,current_k1->data);
                if(mpz_cmp(mod_result1,G_0)==0 && mpz_cmp(mod_result2,G_1)==0){
                    count++; //20 - > 5 and 2 
                }
                current_div=current_div->next;
                mpz_clear(G_d);
                mpz_clear(G_d_1);
                mpz_clear(temp_res);
                mpz_clear(mod_result1);
                mpz_clear(mod_result2);
            }
            if(count==0){
               // gmp_printf("\nthe current k1 mod is %Zd\n", current_k1->data);
                insertAtBeginning(&final_list,current_k1->data);
            }
        }
        current_k1=current_k1->next;
        mpz_clear(result_gcd);
    }
    while(current_otherModuli!=NULL){
        mpz_t gcd_result;       //free this later- i did 
        mpz_init(gcd_result);   
        gcd(det,current_otherModuli->data,gcd_result);
        if(mpz_cmp_ui(gcd_result,1)==0){
            current_otherModuli=current_otherModuli->next;
        }
        else{
            mpz_t mod_result1, mod_result2;    //free those later- i did 
            mpz_init(mod_result1);
            mpz_init(mod_result2);
            mpz_mod(mod_result1,fib_num,current_otherModuli->data);
            mpz_mod(mod_result2,second_fib_num,current_otherModuli->data);
            if(mpz_cmp(mod_result1,G_0)!=0 || mpz_cmp(mod_result2,G_1)!=0){
                current_otherModuli=current_otherModuli->next;
            }
            else{
                int count2=0;
                current_div=divisor;
                while(current_div!=NULL && count2==0){
                    mpz_t fib3, fib4, temp_result, mod_result3, mod_result4;  //free those later- i did 
                    mpz_init(fib3);
                    mpz_init(fib4);
                    mpz_init(temp_result);
                    mpz_init(mod_result3);
                    mpz_init(mod_result4);
                    mpz_divexact(temp_result,converted_period,current_div->data);
                    unsigned long long data= mpz_get_ui(temp_result);
                    generalizedNthFibNum(data,G_0,G_1,fib3);
                    generalizedNthFibNum(data+1, G_0,G_1,fib4);
                    mpz_mod(mod_result3,fib3,current_otherModuli->data);
                    mpz_mod(mod_result4, fib4, current_otherModuli->data);
                    if(mpz_cmp(mod_result3,G_0)==0 && mpz_cmp(mod_result4,G_1)==0){
                        count2++;
                    }
                    current_div=current_div->next;
                    mpz_clear(fib3);
                    mpz_clear(fib4);
                    mpz_clear(temp_result);
                    mpz_clear(mod_result3);
                    mpz_clear(mod_result4);
                }
                if(count2==0){
                    insertAtBeginning(&final_list,current_otherModuli->data);
                }
                current_otherModuli=current_otherModuli->next;
            }
            mpz_clear(mod_result1);
            mpz_clear(mod_result2);
        }
        mpz_clear(gcd_result);
    }
    
    deleteList(&divisor);
    return final_list;
}
//testing if this loop runs
int main() {
    
   /* unsigned long long k_1 = 50;
    struct Node* m1_head = NULL;
    struct Node* div1_head = NULL;
    struct Node* f1_head = NULL;

    mpz_t fib_num;
    mpz_init(fib_num);
    //theNthLucasNum(k_1 , 0, 1, fib_num);
    mpz_t a, b;
    mpz_t c,d;
    mpz_init_set_ui(c,50);
    mpz_init_set_ui(d,29);
    mpz_t result;
    mpz_init(result);
    mpz_init(a);
    mpz_init(b);
    mpz_set_ui(a,0);
    mpz_set_ui(b,1);
   
    generalizedNthFibNum(k_1,a,b,fib_num);
    findDivisors(fib_num, &m1_head);
    displayList(m1_head);
    gmp_printf("\nthe %llu fib num is %Zd\n", k_1,fib_num);
    mpz_clear(fib_num);
    gcd(c,d,result);
    gmp_printf("%Zd is the result\n", result);
    // Further logic as per the original code, adapted for GMP
    mpz_clear(result);
    deleteList(&m1_head);
    deleteList(&f1_head);
    deleteList(&div1_head);

    mpz_t n;
    mpz_init_set_ui(n,101);
    int res= mpz_isPrime(n);
    printf("result is %d\n", res);*/

   // mpz_t period;
    //mpz_init_set_ui(period,6);
    mpz_t G_0, G_1,det, first, second, gcd_result;
    mpz_init(G_0);
    mpz_init(G_1);
    mpz_init(gcd_result);
    mpz_init(first);
    mpz_init(second);
    unsigned long long a,b, period;
    printf("enter the first number: ");
    scanf("%llu",&a);
    printf("enter the second number: ");
    scanf("%llu", &b);
    mpz_set_ui(G_0,a);
    mpz_set_ui(G_1,b);
    unsigned long long converted_det= pow(b,2) - a*b - pow(a,2);
    mpz_init_set_ui(det,converted_det);
    printf("enter the period h(m): ");
    scanf("%llu", &period);
    generalizedNthFibNum(period,G_0,G_1,first);
    mpz_sub(first,first,G_0);
    generalizedNthFibNum(period+1,G_0,G_1,second);
    mpz_sub(second,second,G_1);
    gcd(first,second,gcd_result);
   // gmp_printf("first is %Zd\n second is %Zd\n", first,second);
    if(isPrime(first)==1 && isPrime(second)==1 && mpz_cmp(first,second)!=0){
        printf("moduli do not exist !");
      //  printf("1");
    }
    else if((isPrime(first)==1 || isPrime(second)==1) && (mpz_divisible_p(first,second)==0 || mpz_divisible_p(second,first)==0)){
        printf("moduli do not exist!");
      //  printf("2");
    }
    else if(mpz_cmp_ui(gcd_result,2)<=0 && period >=3){
        printf("moduli do not exist!");
      //  printf("3");
    }
    else{
        if(mpz_cmp_ui(det,1)==0 || mpz_cmp_ui(det,-1)==0){
            if(period%2!=0){
                printf("there is no m's!");
            }
            else{
                struct Node* k1_list= k_1Moduli(period);
                printf("moduli: ");
                displayList(k1_list);
            }
        }
        else{
            struct Node* newNode= potentialModuli(period,converted_det);
            struct Node* k1_list= k_1Moduli(period);
            printf("\nk1 list: ");
            displayList(k1_list);
            printf("\nother moduli: ");
            displayList(newNode);
            struct Node* final_moduli_res= CalculateGenFibMod(k1_list,newNode,G_0,G_1,period,det);
            printf("\nfinal result: ");
            displayList(final_moduli_res);
        }
    }
    mpz_clear(G_0);
    mpz_clear(G_1);
    mpz_clear(det);
    mpz_clear(gcd_result);
    mpz_clear(first);
    mpz_clear(second);
    return 0;
}
