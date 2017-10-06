/*
 * Introduction to HPC 
 * Assignment 3
 * Program: 2
 * Hitender Prakash (hprakash@iu.edu) 
 */
#include <stdio.h>
#include "mpi.h"

typedef struct{
	int max_iter;
	double t0;
	double tf;
	double xmax[12];
	double xmin;
}Pars;

int main(int argc, char *argv[]){
	int myid, numprocs, left, right;
	Pars buffer, buffer2;
	MPI_Request request;
	MPI_Status status;
	
	MPI_Init(&argc, &argv);
	//section for creating derived datatype in MPI - MPI_Pars
	MPI_Datatype MPI_Pars;
	int blocks[]={1,1,1,12,1};
	MPI_Datatype types[]={MPI_INT,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
	MPI_Aint intex,doublex;
	MPI_Type_extent(MPI_INT,&intex);
	MPI_Type_extent(MPI_DOUBLE,&doublex);
	MPI_Aint disp[5];
	disp[0]=0;
	disp[1]=intex;
	disp[2]=intex+doublex;
	disp[3]=intex+doublex+doublex;
	disp[4]=intex+doublex+doublex+doublex;
	MPI_Type_struct(5,blocks,disp,types,&MPI_Pars);
	MPI_Type_commit(&MPI_Pars);
	//MPI_Pars created
	
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	
	right=(myid+2)%numprocs;
	left=myid-2;
	if(left<0){
		left=left+numprocs;
	}
	
	buffer.max_iter=myid;
	buffer.t0=3.14*myid;
	buffer.tf=1.67*myid;
	buffer.xmin=2.55*myid;
	for(int i=0;i<12;i++){
		buffer.xmax[i]=2.7*myid;
	}
	
	MPI_Sendrecv(&buffer,1,MPI_Pars,right,0, &buffer2, 1, MPI_Pars, left, 0, MPI_COMM_WORLD,&status);
	printf("\nProcess %d received %d \n",myid,buffer2.max_iter);
	
	//free the newly created datatype MPI_Pars
	MPI_Type_free(&MPI_Pars);
	MPI_Finalize();
	return 0;
}
