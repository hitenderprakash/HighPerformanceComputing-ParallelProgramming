/* Introduction to High performance Computing 
 * Assignment 2
 * Question 4: Parallel Matrix-Vector Muliplicatin and Strong scaling
 * 
 * Author Hitender Prakash (hprakash@iu.edu)
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>

int main(){
	const int size=10000;
	int i,j;
	
	//Allocating large size arrays can cause program to crash 
	//therefor changing the implementation to use dynamically allocated memory for array
	
	//double A[size*size]; 
	//double x[size];
	//double b[size];
		
	double **A=0,*x=0, *b=0;
	
	//Matrix A can be a one-D array like:  double *A = (double*)malloc((size*size)*sizeof(double));
    //but I like to access matrices like 2D arrays
	A=(double **)malloc(size*sizeof(double *));
	for(i=0;i<size;i++){
		A[i]=(double *)malloc(size*sizeof(double));
	}	
	
	x=(double *)malloc(size*sizeof(double));
	b=(double *)malloc(size*sizeof(double));
	
	for(j=0;j<size;j++){
		for(i=0;i<size;i++){
			//A[i+size*j]=sin(0.01*(i+size*j)); //if one-D array is used for A
			A[j][i]=sin(0.01*(i+size*j));
		}
		b[j]=cos(0.01*j);
		x[j]=0.0;
	}
	
	//int threads; //to store number of threads()
	clock_t start, finish;
	
	start=clock();
	//comment this line for sequential execution
	#pragma omp parallel for shared(j) private(i)
	for(j=0;j<size;j++){
		//threads=omp_get_num_threads();
		for(i=0;i<size;i++){
			//x[j]+= A[i+size*j]*b[i]; //if one-D array is used for A
			x[j]+= A[j][i]*b[i];
		}
	}
	finish=clock();
	
	double time_taken=((double)(finish-start)/CLOCKS_PER_SEC);
	printf("\nx[%d] = %g",5050,x[5050]);
	printf("\nTime taken: %f\n",time_taken);
	
	//free memory block
	for(i=0;i<size;i++){
		free(A[i]);
	}
	free(A);
	free(x);
	free(b);
	//free memory block ends
	return 0;
}
