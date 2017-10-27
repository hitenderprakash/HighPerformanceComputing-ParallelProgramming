/*
 * Introduction to HPC 
 * Program: Calculate PI using MonteCarlo Simulation
 * Hitender Prakash (hprakash@iu.edu) 
 * version 2 
 * Date Oct 25, 2017
 */
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

//to find min
#define min(a,b) (((a)<(b))?(a):(b))

double monteCarloSimulation(long totalSamples){
	//we run the simulation totalSamples time in the first quadrent of unit circle
	//record the ration of points which falls in circle 
	//to the total points
	double pi_from_one_simulation=0.0;
	long count=0; //count of points which falls inside the circle
	long i=0; //loop counter
	for(i=0;i<totalSamples;i++){
		//generate x cordinate
		double x=(double)rand()/(double)RAND_MAX; //keep it between 0 and 1
		//generate random y coordinate
		double y=(double)rand()/(double)RAND_MAX; //keep it between 0 and 1
		//see if it falls in the unit circle area
		if((x*x+y*y)<=1){
			count++;
		}
	}
	pi_from_one_simulation=4.0*(double)count/(double)totalSamples;
	return pi_from_one_simulation;
}

int main(int argc,char **argv) {
  //pulling this code up so that each process share the same vectors 
  // Make the local vector size constant
  int total_number_of_simulations = 10000;
  //initialize the vectors
       
  MPI_Init(&argc,&argv);
  int rank,p,i, root = 0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&p);


  //section to divide the task between processes
  int chunk=0; //chunk size to be taken up by single process
  if(total_number_of_simulations%p==0){
	  chunk=(total_number_of_simulations/p);
  }
  else{
	  chunk=(total_number_of_simulations/p)+1;
  }
  //decide loop boundary to be executed by process 
  int startLoop=chunk*rank;
  int endLoop=min(startLoop+chunk-1,total_number_of_simulations-1);

  double localSum=0.0;
  for (i=startLoop;i<=endLoop;i++) {	  
	  double Pi_from_one_simulation=monteCarloSimulation(2000000000);
	  localSum += Pi_from_one_simulation;    
  } 
  //reduction to sum them all
  double Sum_of_pi=0.0;
  MPI_Reduce(&localSum,&Sum_of_pi,1,MPI_DOUBLE,MPI_SUM,root,MPI_COMM_WORLD);
  if ( rank == root ) {
    printf("\nThe Final value of PI: %.10lf",Sum_of_pi/(double)total_number_of_simulations);
  }
  
  MPI_Finalize();
  return 0;
}


