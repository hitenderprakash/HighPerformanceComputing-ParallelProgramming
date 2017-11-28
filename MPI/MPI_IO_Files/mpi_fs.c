/*
 * Program to demonstarte MPI Paralle I/O
 * create, write and read parallel file
 * Hitender Prakash
 * 
 */

#include <stdio.h>
#include <mpi.h>

int main(int argc, char **argv){
	int i;
	int rank;
	int size;
	int offset;
	int nints;
	int N=20;
	
	MPI_File fhw;
	MPI_Status status;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	int buf[N];
	for(i=0;i<N;i++){
		buf[i]=rank*N+i;
	}
	
	MPI_File_open(MPI_COMM_WORLD, "datafile", MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fhw);
	
	//each process write the "buf" array to the file called "datafile"
	offset=rank*N*sizeof(int);
	MPI_File_write_at(fhw,offset,buf,N,MPI_INT,&status);	
	//close the file
	MPI_File_close(&fhw);
	
	//reopen the file
	MPI_File_open(MPI_COMM_WORLD, "datafile", MPI_MODE_RDONLY,MPI_INFO_NULL,&fhw);
	//read the integer array back and print first four elements
	MPI_File_read_at(fhw,offset,buf,N,MPI_INT,&status);
	printf("Process %d read %d %d %d %d\n",rank,buf[0],buf[1],buf[2],buf[3]);
	MPI_File_close(&fhw);
	
	MPI_Finalize();
	return 0;
}
