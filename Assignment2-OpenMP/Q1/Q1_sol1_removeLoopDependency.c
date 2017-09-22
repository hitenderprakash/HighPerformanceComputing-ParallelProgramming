/*
 * Introduction to High Performance Computing 
 * Assignment 2
 * Question 1
 * 
 * The problem with code was that "dot_prod" variable is not thread safe
 * Solution 1: To remove loop dependency instead of using a single shared variable take an array  so that each
 * thread write to separate array element without any race condition
 * Later add all the array element and that will be the dot_product
 */

#include <stdio.h>
#include <omp.h>

int main(){
	const int N=100;
	int i,k;
	float a[N],b[N];	
	float dot_prod[N];
	
	for(k=0;k<N;k++){
		a[k]=3.0*k;
		b[k]=1.8*k;
	}
	
	#pragma omp parallel 
	{
		#pragma omp for
		for(i=0;i<N;i++){
			dot_prod[i]=a[i]*b[i]; //each thread write its corresponding element without synchronization issue
		}
	}
	
	float final_dot_prod=0.0; // to store the sum of all array elemnt 
	for(i=0;i<N;i++){final_dot_prod+=dot_prod[i];}
	printf("\nInner product of a[] and b[] = %f\n", final_dot_prod);
	return 0;
}

