/*
 * Introduction to High Performance Computing 
 * Assignment 2
 * Question 2 - Scheduling - static, dynamic
 */
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>

int main(){
	const int N=28;
	int nthreads, threadid,i;
	double a[N],b[N],result[N];
	
	//initialize
	for(i=0;i<N;i++){
		a[i]=(1.0)*(float)i;
		b[i]=(2.0)*(float)i;
	}

	int chunk=7;
	#pragma omp parallel private(threadid,i)
	{
		threadid=omp_get_thread_num();
		
		//for static scheduling
		#pragma omp for schedule(static,chunk)
		
		//for dynamic scheduling
		//#pragma omp for schedule(dynamic)
		for(i=0;i<N;i++){
			result[i]=a[i]+b[i];
			printf("\nThread id: %d working on index: %d", threadid,i);
		}
	}
	printf("\nTest Results[19] = %g\n", result[19]);	
	return 0;
}
