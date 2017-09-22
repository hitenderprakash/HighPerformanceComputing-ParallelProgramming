/*
 * Introduction to High Performance Computing 
 * Assignment 2
 * Question 3 - Sections
 */
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>

int main(){
	int N=1000;
	int x[N],i=0,max_x=0,min_x=100,sum=0,sum2=0;
	float mean=0.0,mean2=0.0;
	long double var=0.0; // to save from overflow

	/*initialize*/
	srand(1.0);//initialize random variable seed
	#pragma omp parallel for
	for(i=0;i<N;i++){
		x[i]=rand();
	}
	
	#pragma omp parallel private(i) shared (x)
	{
		#pragma omp sections
		{
			/*fork 3 different threads */
			#pragma omp section
			{
				//printf("\nThread: %d started Section 1",omp_get_thread_num());
				for(i=0;i<N;i++){
					if(x[i]>max_x){max_x=x[i];}
					if(x[i]<min_x){min_x=x[i];}
				}
				printf("\nThe max of x = %d", max_x);
				printf("\nThe min of x = %d", min_x);
				//printf("\nThread: %d finished Section 1",omp_get_thread_num());
			}
			#pragma omp section
			{
				//printf("\nThread: %d started Section 2",omp_get_thread_num());
				for(i=0;i<N;i++){
					sum=sum+x[i];
				}
				mean=sum/N;
				printf("\nMean of x = %f", mean);
				//printf("\nThread: %d finished Section 2",omp_get_thread_num());
			}
			#pragma omp section
			{
				//printf("\nThread: %d started Section 3",omp_get_thread_num());
				for(i=0;i<N;i++){
					sum2=sum2+x[i]*x[i];
				}
				mean2=sum2/N;
				printf("\nMean2 of x = %f", mean2);
				//printf("\nThread: %d finished Section 3",omp_get_thread_num());
			}
		}
	}
	var=(long double)mean2-(long double)(mean*mean); //type cast to long double otherwise value will overflow
	printf("\nVarience of x = %llf\n", var);	
	return 0;
}
