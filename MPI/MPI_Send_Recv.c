#include <stdio.h>
#include "mpi.h"
#include <string.h>

int main(int argc, char *argv[]){
	int myrank,size,source,dest,tag=0;
	char msg[100];
	MPI_Status status;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	
	if(myrank !=0){
		dest=0;
		printf("\nProcess: %d sending message to Process: 0",myrank);
		sprintf(msg,"[Hello %d]",myrank);
		MPI_Send(msg,strlen(msg)+1,MPI_CHAR,dest,tag,MPI_COMM_WORLD);
	}
	else{
		for(source=1;source<size;source++){
			MPI_Recv(msg,100,MPI_CHAR,source,tag,MPI_COMM_WORLD,&status);
			printf("\nProcess: 0 receiving message %s from Process:%d",msg,source);
		}
	}
	MPI_Finalize();
	printf("\n");
	return 0;
}
