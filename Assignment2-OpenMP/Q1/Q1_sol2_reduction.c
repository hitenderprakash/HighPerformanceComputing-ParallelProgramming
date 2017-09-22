/*
 * Introduction to High Performance Computing 
 * Assignment 2
 * Question 1
 * The problem with code was that "dot_prod" variable is not thread safe
 * Solution 2: Use REDUCTION
 */

#include <stdio.h>
#include <omp.h>

int main(){
	const int N=100;
	int i,k;
	float a[N],b[N];	
	float dot_prod;
	dot_prod=0.0;
	
	for(k=0;k<N;k++){
		a[k]=3.0*k;
		b[k]=1.8*k;
	}
	//reduction with + operator i.e for each thread compiler makes a temporary copy of dot_prod initialized with identity value for + i.e 1
	//And in the end it adds all the value 
	#pragma omp parallel for reduction(+:dot_prod) 
		for(i=0;i<N;i++){
			dot_prod+=a[i]*b[i]; 
		}
	printf("\nInner product of a[] and b[] = %f\n", dot_prod);
	return 0;
}

