/*
 * Hitender Prakash 
 * Program: Jacobi iterative method (Heat transfer simulation)
 * Parallel implementation using MPI
 * 
 * This version of program takes lesser memory as the original matrix is maintained by proc 0 only
 * rest of the procs keep only the smaller chunk of the matrix.
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "mpi.h"

#define COLS 21
#define ROWS 21
#define TEMP 10.0 //Temperature at source cell
#define DEBUG 1   //debug flag, shows intermediate transitions
#define EPS 1e-6  //minimum change in temp. between trasitions for quitting
#define I_FIX 9   //Temperature source coordinates
#define J_FIX 9
#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

//Method for printing the sheet (matrix)
void print_matrix(double** matrix,int rows, int cols){
    for (int i = 0; i < rows; i++) {
		printf("\n");
        for (int j = 0; j < cols; j++){
            printf("%.6f	", matrix[i][j]);
        }
    }
}

//prints matrix to console and write to a file "results.txt" in current directory
void print_and_write_matrix(double** matrix,int rows, int cols){
	FILE *fout=fopen("results.txt","a");
    for (int i = 0; i < rows; i++) {
		printf("\n");
		fprintf(fout,"\n");
        for (int j = 0; j < cols; j++){
            printf("%.6f	", matrix[i][j]);
            fprintf(fout,"%.6f	", matrix[i][j]);
        }
    }
    fclose(fout);
}

//init matrix with default value as zero
double** alloc_matrix(int rows, int cols){
    double** matrix = (double**) malloc(rows * sizeof(double *));
    int i;
    for(i=0;i<rows;i++){matrix[i] = (double*) calloc(cols,sizeof(double));} 
    return matrix;
}

//Method for computing the temperature of a cell from its surrounding cells
void computeJacobi(double **oldmat, double **newmat, int chunk, int cols, double *ghostrow_recv_upper, double *ghostrow_recv_lower){	
	int i,j;
	for(i=0;i<chunk;i++){
		for(j=0;j<cols;j++){
			double upper=0.0;
			double lower=0.0;
			double right=0.0;
			double left=0.0;
			
			if(j!=0){left=oldmat[i][j-1];}
			if(j!=cols-1){right=oldmat[i][j+1];}
			
			if(i!=0){upper=oldmat[i-1][j];}			
			else if(i==0){upper=ghostrow_recv_upper[j];} 
			
			if(i!=chunk-1){lower=oldmat[i+1][j];}
			else if(i==chunk-1){lower=ghostrow_recv_lower[j];}		
			newmat[i][j]= (0.25*(upper+lower+right+left));
		}
	}
}

//copying two equal size matrices
void copyMatrix(double **dest, double **source,int rows, int cols){
	int i,j;
	for(i=0;i<rows;i++){
		for(j=0;j<cols;j++){
			dest[i][j]=source[i][j];
		}
	}
}

//find the max temp. difference between any two cells after one transition
double findMaxDiff(double **oldMat, double **newMat, int chunk, int cols){
	double max_diff=DBL_MIN;
	int i,j;
	for(i=0;i<chunk;i++){
		for(j=0;j<cols;j++){
			max_diff= max(max_diff,fabs(newMat[i][j]-oldMat[i][j]));
		}
	}
	return max_diff;
}

// Main program
int main(int argc, char *argv[]) {
	
	MPI_Status status;
	int rows=ROWS;
	int cols=COLS;
	
	double **a_old;
	double **subMat;
	double **oldsubMat;
	
	double maxerr;
	double global_err;
	double old_global_err;
	
	double *ghostrow_send_upper=(double*) calloc(cols,sizeof(double));
	double *ghostrow_send_lower=(double*) calloc(cols,sizeof(double));
	double *ghostrow_recv_upper=(double*) calloc(cols,sizeof(double));
	double *ghostrow_recv_lower=(double*) calloc(cols,sizeof(double));
		
    MPI_Init(&argc,&argv);
	int rank,numProc, root = 0;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&numProc);
	//printf("\nRank: %d, Number of Processes: %d.",rank,numProc);
    
    int chunk=0; //rows to be taken up by single process
	if(rows%numProc==0){
		chunk=(rows/numProc);
	}
	else{
		chunk=(rows/numProc)+1;
	}
	//printf("\nEach proc gets: %d rows",chunk);
	int startRow=rank*chunk;
	
	int paddedEndRow=rank*chunk+chunk-1;
	int endRow=min(paddedEndRow,rows-1); 
	//end row will not contain any padded row for example 21 rows matrix was given to 4 procs
	//each proc will have 6 rows i.e 21 row's matrix will be padded with 3 more rows i.e augmented matrix has 24 rows
	//but for each proc endRow will have valid row from the actual matrix i.e 4th proc will have end row - 21 though it can also access padded rows as well
	//normally every process except last one will have same value for paddedEndRow and endRow or if no padding was required then last proc will also have same value for both
	
	//printf("\nRank %d gets from %d to %d", rank, startRow, endRow);
	
	//for mapping the temp. source to proc when rows gets divided between different procs
	int proc_at_Fixed_Temp=(I_FIX/chunk);
	int Fixed_row=I_FIX%chunk;

	MPI_Barrier(MPI_COMM_WORLD);
	//printf("\nMy Rank: %d, Fixed Proc: %d  Fixed Row: %d",rank,proc_at_Fixed_Temp,Fixed_row);
	
	//now need to alloc and init the matrix with perfect size (padding may be required) 
	//this will be done only for proc 0
	if(rank==0){
		a_old=alloc_matrix(chunk*numProc, cols);
		//set the temperature at source
		a_old[I_FIX][J_FIX] = TEMP;
	}
	
	
	//initialize sub matrix with less number of arrays
	//in this matrix each process will keep newly computed temperature values
	//This matrix is small in size and will keep data for one process only
	oldsubMat = alloc_matrix(chunk, cols); //proc 0 divide matrix in this sub matrix
	subMat = alloc_matrix(chunk, cols);    //used for newly computed values
	
	//proc 0 distribute the matrix between all other processes
	if(rank==0){
		int i,j;
		//send to other procs
		for(i=1;i<numProc;i++){
			for(j=0;j<chunk;j++){
				MPI_Send(a_old[i*chunk+j], cols, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			}
		}
		//copy data to its own oldsubMat
		for(i=0;i<chunk;i++){
			for(j=0;j<cols;j++){
				oldsubMat[i][j]=a_old[i][j];
			}
		}
	}
	else{
		int i,j;
		for(j=0;j<chunk;j++){
			MPI_Recv(oldsubMat[j], cols, MPI_DOUBLE, 0,0 , MPI_COMM_WORLD, &status);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	
	while(1){
		
		if(rank==0){
			int i=0;
			for(i=0;i<cols;i++){ghostrow_send_lower[i]=oldsubMat[chunk-1][i];}
			// will have to take care here that if there is only one processor available there will be no next proc
			//otherwise send and receive will fail
			if(numProc>1){
				MPI_Send(ghostrow_send_lower, cols, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
				MPI_Recv(ghostrow_recv_lower, cols, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &status);
			}
		}
		if(rank==numProc-1){
			int i=0;
			for(i=0;i<cols;i++){ghostrow_send_upper[i]=oldsubMat[0][i];}
			// will have to take care here that if there is only one processor available there will be no previous proc
			//otherwise send and receive will fail
			if(numProc>1){
				MPI_Send(ghostrow_send_upper, cols, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
				MPI_Recv(ghostrow_recv_upper, cols, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
			}
		}
		if(rank!=0 && rank!=numProc-1){
			int i=0;
			for(i=0;i<cols;i++){ghostrow_send_upper[i]=oldsubMat[0][i];}
			for(i=0;i<cols;i++){ghostrow_send_lower[i]=oldsubMat[chunk-1][i];}
			
			MPI_Send(ghostrow_send_lower, cols, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
			MPI_Recv(ghostrow_recv_lower, cols, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &status);
			MPI_Send(ghostrow_send_upper, cols, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
			MPI_Recv(ghostrow_recv_upper, cols, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
		}		

		//Computer new temperature values from the surrounding cells
		computeJacobi(oldsubMat, subMat, chunk,cols,ghostrow_recv_upper,ghostrow_recv_lower);
		
		
		//======================only required if all side cells to be kept 0==================
		//this loop will make all element of first and last coloumn zero in old large matrix 
		//as well as in its mapping small matrix at particular proc (run for all proc)
		
		//it could have been achieved in ComputeJacobi method by altering inner loop from j=1 to j=cols-2 instead of j=0 to j=cols-1
		//but wanted to keep that code general, so doing it separately here
		int i,j;
		for(i=0;i<chunk;i++){
			subMat[i][0]=0.0;
			subMat[i][cols-1]=0.0;
			oldsubMat[i][0]=0.0;
			oldsubMat[i][cols-1]=0.0;
		}
		//this loop will make all element of first row zero in old large matrix 
		//as well as in its mapping small matrix at particular proc (always proc 0)
		if(rank==0){
			for(i=0;i<cols;i++){
				subMat[startRow][i]=0.0;
				oldsubMat[startRow][i]=0.0;
			}
		}
		//this loop will make all element of last row zero in old large matrix 
		//as well as in its mapping small matrix at particular proc (normally last proc)
		//as well as make all padded rows elemnt zero 
		if(paddedEndRow >= rows-1){
			for(i=paddedEndRow;i>=startRow && i>=rows-1;i--){ //also need to make sure that we anly loop till startRow, after that rows belong to another procs
				for(j=0;j<cols;j++){
					subMat[i-startRow][j]=0.0; 
					oldsubMat[i-startRow][j]=0.0; 
				}
			}
		}
		//======================only required if all side cells to be kept 0==================
			
		//zero out the extra padded rows in original large matrix.
		//also need to do it in coreesponding small matrix taken for computing new values so that ...
		//..does not impact the computation of maxerr
		/*int i,j;				
		if(paddedEndRow > rows-1){
			for(i=paddedEndRow;i>=startRow && i > rows-1;i--){ //also need to make sure that we anly loop till startRow, after that rows belong to another procs
				for(j=0;j<cols;j++){
					subMat[i-startRow][j]=0.0; 
					oldsubMat[i-startRow][j]=0.0; 
				}
			}
		}
		*/
		//adjust the fixed point in new small matrix
		if(rank==proc_at_Fixed_Temp){
			subMat[Fixed_row][J_FIX]=TEMP;
			oldsubMat[Fixed_row][J_FIX]=TEMP;
		}
		
		//compute the max diff here for every thread
		maxerr=findMaxDiff(oldsubMat,subMat,chunk,cols);
		//printf("\n The %d 's Maxerr is: %lf",rank,maxerr);
		
		//copy back the subMatCP to a_old corresponding elements 		
		copyMatrix(oldsubMat,subMat,chunk,cols);
		
		
		//AT proc 0 reconstruct the entire matrix by taking new values from each process
		/*
		if(rank==0){
			int i;
			for(i=1;i<numProc;i++){MPI_Recv(a_old[i*chunk], cols*chunk, MPI_DOUBLE, i,0 , MPI_COMM_WORLD, &status);}
		}
		else{MPI_Send(a_old[rank*chunk], cols*chunk, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);}
		*/
		
		if(rank==0){
			int i,j;
			for(i=1;i<numProc;i++){
				for(j=0;j<chunk;j++){
					MPI_Recv(a_old[i*chunk+j], cols, MPI_DOUBLE, i,0 , MPI_COMM_WORLD, &status);
				}
			}
			//copy data back from its own oldsubMat
			for(i=0;i<chunk;i++){
				for(j=0;j<cols;j++){
					a_old[i][j]=oldsubMat[i][j];
				}
			}
		}
		else{
			int i;
			for(j=0;j<chunk;j++){
				MPI_Send(oldsubMat[j], cols, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		
		//print if debugging is enabled
		if(DEBUG){
			if(rank==0){
				printf("\n\n");
				print_matrix(a_old,rows,cols);
			}

		}		
		MPI_Barrier(MPI_COMM_WORLD);
		
		//find the max error among all processes
		MPI_Reduce(&maxerr,&global_err,1,MPI_DOUBLE,MPI_MAX,root,MPI_COMM_WORLD);
		MPI_Bcast(&global_err, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		//printf("\n The global_err is: %lf",global_err);

		if(global_err<EPS){
			break;
		}
		//In case error does not change is stays fixed but larger than EPS
		if(global_err==old_global_err){break;}
		old_global_err=global_err;
	}
	//print final state of matrix (sheet)
	if(rank==0){
		//final display
		printf("\n\nFinal Matrix:\n");
		//print_matrix(a_old,numProc*chunk,cols); //print padded rows also 
		//print_matrix(a_old,rows,cols);  //do not print padded rows
		
		//print matrix and also write it to file "results.txt"
		print_and_write_matrix(a_old,rows,cols);
		printf("\n\n");
	}
	//FREE MEMORY
	int i;
	for(i=0;i<chunk;i++){free(subMat[i]);}
	free(subMat);
	for(i=0;i<chunk;i++){free(oldsubMat[i]);}
	free(oldsubMat);
	if(rank==0){
		for(i=0;i<chunk*numProc;i++){free(a_old[i]);}
		free(a_old);
	}
	free(ghostrow_send_upper);
	free(ghostrow_send_lower);
	free(ghostrow_recv_upper);
	free(ghostrow_recv_lower);
	//FREE MEMORY
	MPI_Finalize();
    return 0;
}
