/* Program to integrate cos(x)*sin(x/2) from 0 to Pi/2
 * Hitender Prakash (hprakash@iu.edu)
 * Version 1: Parallel code with MPI
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h> //link with -lm compiler flag
#include "mpi.h"

#define PI 3.14159265

#define min(a,b) (((a)<(b))?(a):(b))
int main(int argc,char **argv) {
  
  double c=0.0;
  double incr=0.00001; //step increment
  double upperLimit=PI/2;
  int reps=(PI/2)/incr +1;
  
  //debug info
  printf("\nTotal reps will be: %d", reps);
  int count=0;
  double integral=0.0;
	    
  MPI_Init(&argc,&argv);
  int rank,p,i, root = 0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&p);


  //section to divide the task between processes
  int chunk=0; //chunk size to be taken up by single process
  if(reps%p==0){
	  chunk=(reps/p);
  }
  else{
	  chunk=(reps/p)+1;
  }
  
  //debug info 
  printf("\nChunk Size: %d",chunk);
  
  
  //decide loop boundary to be executed by process 
  double startLoop=chunk*rank*incr;
  double endLoop=min(startLoop+chunk*incr,upperLimit); //make sure do not overflow upper limit
  double theta=startLoop;
  
  double localSum=0.0;
  
  //debug info 
  printf("\Rank: %d- I will do from: %lf to %lf", rank, startLoop, endLoop);
  
  while(theta <=endLoop){
	  localSum+=cos(theta)*sin(theta/2);
	  theta+=incr;
  }
  //reduction to sum them all
  double Sum=0.0;
  MPI_Reduce(&localSum,&Sum,1,MPI_DOUBLE,MPI_SUM,root,MPI_COMM_WORLD);
  if ( rank == root ) {
    printf("\nThe Final value of integration: %.10lf",Sum*incr);
  }
  
  MPI_Finalize();
  return 0;
}


