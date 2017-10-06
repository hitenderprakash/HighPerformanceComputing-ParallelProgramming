#include <stdio.h>
#include "mpi.h"

int main(int argc, char *argv[]){
	
	//initialize MPI constructs
	MPI_Init(&argc, &argv);
	printf("\nWelcome to the world of MPI");
	MPI_Finalize();
	return 0;
}
