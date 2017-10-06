/*
 * Introduction to HPC 
 * Assignment 3
 * Program: 3
 * Hitender Prakash (hprakash@iu.edu) 
 */
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h> //add flag -lm while compiling

//to find min
#define min(a,b) (((a)<(b))?(a):(b))

int main(int argc,char **argv) {
  //pulling this code up so that each process share the same vectors 
  // Make the local vector size constant
  int global_vector_size = 10000;
  // initialize the vectors
  double *a, *b;
  a = (double *) malloc(global_vector_size*sizeof(double));
  b = (double *) malloc(global_vector_size*sizeof(double));
       
  MPI_Init(&argc,&argv);
  int rank,p,i, root = 0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&p);

  double pi = 4.0*atan(1.0);

  //section to divide the task between processes
  int chunk=0; //chunk size to be taken up by single process
  if(global_vector_size%p==0){
	  chunk=(global_vector_size/p);
  }
  else{
	  chunk=(global_vector_size/p)+1;
  }
  //decide loop boundary to be executed by process 
  int startLoop=chunk*rank;
  int endLoop=min(startLoop+chunk-1,global_vector_size-1);
  
  for (i=startLoop;i<=endLoop;i++) {
    a[i] = sqrt(i); 
    b[i] = sqrt(i); 
  }
  // compute the dot product
  double localSum = 0.0;
  for (i=startLoop;i<=endLoop;i++) {
    localSum += a[i]*b[i];
  }
  
  //reduction to sum them all
  double Sum=0.0;
  MPI_Reduce(&localSum,&Sum,1,MPI_DOUBLE,MPI_SUM,root,MPI_COMM_WORLD);
  if ( rank == root ) {
    printf("The dot product is %g.  Answer should be: %g\n",Sum,0.5*global_vector_size*(global_vector_size-1));
  }
  
  MPI_Finalize();
  
  free(a);
  free(b);
  return 0;
}


