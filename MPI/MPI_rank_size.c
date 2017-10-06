#include <stdio.h>
#include "mpi.h"

int main(int argc, char *argv[]){
	
	//initialize MPI constructs
	int size,rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	printf("\nHello from process: %d of total %d processes\n",rank,size);
	MPI_Finalize();
	return 0;
}
