/*
 * Introduction to High Performance Computing 
 * Assignment 2
 * Question 1
 * Input buggy program
 * The problem with code is that "dot_prod" variable is not thread safe
 * 
 * version 1
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
	
	#pragma omp parallel 
	{
	#pragma omp for
		for(i=0;i<N;i++){
			dot_prod=dot_prod +a[i]*b[i];
		}
	}

	printf("\nInner product of a[] and b[] = %f\n", dot_prod);
	return 0;
}
