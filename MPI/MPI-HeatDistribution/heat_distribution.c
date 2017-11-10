#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "mpi.h"

#define COLS 12
#define ROWS 12
#define TEMP 50.0
#define DEBUG 1
#define EPS 1e-6
#define I_FIX 7
#define J_FIX 5

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

double max_abs(double** m1, double** m2, int rows, int cols){
    double max_val = DBL_MIN;
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++){
            if (fabs(m1[i][j] - m2[i][j]) > max_val) {
                max_val = fabs(m1[i][j] - m2[i][j]);
            }
        }
    return max_val;
}

void print_matrix(double** matrix,int rows, int cols){
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++)
            printf("%f ", matrix[i][j]);
        printf("\n");
    }
}

void copy_matrix(double** dest, double** source,int rows, int cols) {
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            dest[i][j] = source[i][j];
}

double** alloc_matrix(int rows, int cols){
    double** matrix = (double**) malloc(rows * sizeof(double *));
    int i;
    for(i=0;i<rows;i++){matrix[i] = (double*) calloc(cols,sizeof(double));} 
    return matrix;
}

void compute_new_values(double** old_matrix, double** new_matrix){
    for (int i = 1; i < ROWS-1; i++)
        for (int j= 1; j < COLS-1; j++)
            new_matrix[i][j] =
                    0.25 * (old_matrix[i-1][j] + old_matrix[i+1][j] +
                            old_matrix[i][j-1] + old_matrix[i][j+1]);
    new_matrix[I_FIX][J_FIX] = TEMP;
}

void init_matrix(double** matrix, int rows, int cols){
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = 0.0;
        }
    matrix[I_FIX][J_FIX] = TEMP;
}

//do the computation of avarage value from surrounding values
/*
void computeJacobi(double **oldmat, double **newmat, int rows, int cols, double *ghostrow_recv_upper, double *ghostrow_recv_lower){
	int i,j;
	for(i=0;i<rows;i++){
		for(j=0;j<cols;j++){
			int upper=0;
			int lower=0;
			int right=0;
			int left=0;
			if(j!=0){left=old[i][j-1];}
			if(j!=cols-1){right=old[i][j+1];}
			if(i!=0){upper=old[i-1][j];}
			else if(i==0){upper=ghostrow_recv_upper[j];}
			if(i!=rows-1){lower=old[i+1][j];}
			else if(i==rows-1){upper=ghostrow_recv_lower[j];}			
			newmat[i][j]= (0.25*(upper+lower+right+left));
		}
	}

}*/

void computeJacobi(double **oldmat, double **newmat, int rows, int cols, double *ghostrow_recv_upper, double *ghostrow_recv_lower,int rank,int chunk){
	int startRow=chunk*rank;	
	int i,j;
	for(i=0;i<chunk;i++){
		for(j=0;j<cols;j++){
			int upper=0;
			int lower=0;
			int right=0;
			int left=0;
			
			if(j!=0){left=oldmat[i+startRow][j-1];}
			if(j!=cols-1){right=oldmat[i+startRow][j+1];}
			
			if(i+startRow!=0){upper=oldmat[i+startRow-1][j];}			
			//else if(i+startRow==0){upper=ghostrow_recv_upper[j];} //not required, it will anyway keep upper zero, beacuase for first row for proc 0 ghost cells from above are all zeros
			
			if(i+startRow!=rows-1){lower=oldmat[i+startRow+1][j];}
			//else if(i+startRow==rows-1){upper=ghostrow_recv_lower[j];}// not required, because for last row of proc n-1 ghost cell from below are also zeros			
			newmat[i][j]= (0.25*(upper+lower+right+left));
		}
	}

}

void copyRow(double *dest, double *source, int size){
	int i;
	for(i=0;i<size;i++){
		dest[i]=source[i];
	}

}

void CopyPartOfMatrix(double **dest, double **source, int size,int rank,int chunk){
	//copying a part (few rows) from a bigger source matrix to smaller destination matrix
	int startRow=rank*chunk;
	int i,j;
	for(i=0;i<chunk;i++){
		for(j=0;j<size;j++){
			dest[i][j]=source[i+startRow][j];
		}
	}
	
}
void copyPartToMatrix(double **dest, double **source, int size,int rank,int chunk){
	//copying a smaller source matrix to a bigger destination matrix as its part
	int startRow=rank*chunk;
	int i,j;
	for(i=0;i<chunk;i++){
		for(j=0;j<size;j++){
			dest[i+startRow][j]=source[i][j];
		}
	}
}

void copy2dTo1d(double *dest, double **source,int size, int rank, int chunk){
		int startRow=rank*chunk;
		int endRow=rank*chunk+chunk-1;
		int i,j;
		for(i=0;i<chunk;i++){
			for(j=0;j<size;j++){			
				dest[i*size+j]=source[i+startRow][j];
				//*dest=source[i][j];
				//dest++;
			}
		}
}

void copy1dTo2d(double **dest, double *source,int size, int rank, int chunk){
		int startRow=rank*chunk;
		int endRow=rank*chunk+chunk-1;
		int i,j;
		for(i=0;i<chunk;i++){
			for(j=0;j<size;j++){
				dest[i+startRow][j]=source[i*size+j];
				//dest[i][j]=*source;
				//source++;
			}
		}
}

double findMaxDiff(double **oldMat, double **newMat, int cols, int rank, int chunk){
	int startRow=rank*chunk;
	double max_diff=DBL_MIN;
	int i,j;
	for(i=0;i<chunk;i++){
		for(j=0;j<cols;j++){
			max_diff= max(max_diff,newMat[i][j]-oldMat[i+startRow][j]);
		}
	}
	return max_diff;
}

void print_array(double *arr, int size){
	printf("\n");
	int i=0;
	for(i=0;i<size;i++){
		printf("%lf ",arr[i]);
	}
}

int main(int argc, char *argv[]) {
	
	MPI_Status status;
	int rows=ROWS;
	int cols=COLS;
	
	double **a_old;
	double **subMat;
	double **subMatCP;
	
	double maxerr=DBL_MAX;
	double global_err=DBL_MAX;
	
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
	int endRow=min(rank*chunk+chunk-1,cols-1);
	//printf("\nRank %d gets from %d to %d", rank, startRow, endRow);
		
	//now need to alloc and init the matrix with perfect size (padding may be required)
	a_old=alloc_matrix(chunk*numProc, cols);
	a_old[I_FIX][J_FIX] = TEMP;
	//initially all processes have the same view of the heat sheet the actual averaging starts withing the loop 
	
	//initialize sub matrix with less number of arrays
	subMat=alloc_matrix(chunk, cols);
	subMatCP=alloc_matrix(chunk, cols);
	
	
	//DEBUG
	while(1){
		//working
		MPI_Reduce(&maxerr,&global_err,1,MPI_DOUBLE,MPI_MAX,root,MPI_COMM_WORLD);
		MPI_Bcast(&global_err, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if(global_err<EPS){
			break;
		}
		
		//rank 0 sends the orig matrix to all procs. Actually this part might not be required
		if(rank==0){
			int i;
			for(i=1;i<numProc;i++){
				MPI_Send(a_old[i*chunk], cols*chunk, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			}
		}
		else{
			MPI_Recv(a_old[rank*chunk], cols*chunk, MPI_DOUBLE, 0,0 , MPI_COMM_WORLD, &status);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		// this communication is not required
		
		//DEBUG1////////////////////////////////////////////////////////////////////
		if(rank==0){
			int i=0;
			for(i=0;i<cols;i++){ghostrow_send_lower[i]=a_old[endRow][i];}
			
			MPI_Send(ghostrow_send_lower, cols, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
			MPI_Recv(ghostrow_recv_lower, cols, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &status);
		}
		else if(rank==numProc-1){
			int i=0;
			for(i=0;i<cols;i++){ghostrow_send_upper[i]=a_old[startRow][i];}
			
			MPI_Send(ghostrow_send_upper, cols, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
			MPI_Recv(ghostrow_recv_upper, cols, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
		}
		else{
			int i=0;
			for(i=0;i<cols;i++){ghostrow_send_upper[i]=a_old[startRow][i];}
			for(i=0;i<cols;i++){ghostrow_send_lower[i]=a_old[endRow][i];}

			MPI_Send(ghostrow_send_lower, cols, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
			MPI_Recv(ghostrow_recv_lower, cols, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &status);
			MPI_Send(ghostrow_send_upper, cols, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
			MPI_Recv(ghostrow_recv_upper, cols, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
		}		
		//working till here
		
		//now do the computation
		computeJacobi(a_old, subMat, chunk*numProc,cols,ghostrow_recv_upper,ghostrow_recv_lower,rank,chunk);
		//compute the max diff here for every thread
		maxerr=findMaxDiff(a_old,subMat,cols,rank,chunk);
		
		//copy back the subMatCP to a_old corresponding elements 		
		copyPartToMatrix(a_old,subMat,cols,rank,chunk);
		
		//readjust fixed element
		a_old[I_FIX][J_FIX] = TEMP;
		
		//also can make padded rows all zero
		int i,j;
		for(i=rows;i<chunk*numProc;i++){
			for(j=0;j<cols;j++){
				a_old[i][j]=0.0;
			}
		}
		
		//now we can copy corresponding parts to proc 0's original matrix
		
		if(rank==0){
			int i;
			for(i=1;i<numProc;i++){
				MPI_Recv(a_old[i*chunk], cols*chunk, MPI_DOUBLE, i,0 , MPI_COMM_WORLD, &status);
			}
		}
		else{
			MPI_Send(a_old[rank*chunk], cols*chunk, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if(rank==0){
			printf("\n\n");
			print_matrix(a_old,rows,cols);
		}
		
		MPI_Barrier(MPI_COMM_WORLD);
		
		//DEBUG1/////////////////////////////////////////////////////////////////////
				
		//break;//only for testing remove otherwise
	}
	/*
	if(rank==0){
		//final display
		printf("\nFinal Matrix: \n");
		print_matrix(a_old,rows,cols);
	}	*/
	MPI_Finalize();
    return 0;
	
	//DEBUG	
	//=============================================================
	/*(all procs send to proc 0)
	///send and receive in loop in case the above method to send multiple rows fails
	* 
	if(rank==0){
		int i,j;
		for(i=1;i<numProc;i++){
			startRow=i*chunk;
			for(j=0;j<chunk;j++){
				MPI_Recv(a_old[i*chunk+j], cols, MPI_DOUBLE, i,0 , MPI_COMM_WORLD, &status);
			}
		}
	}
	else{
		int i;
		int startRow=rank*chunk;
		for(j=0;j<chunk;j++){
			MPI_Send(a_old[rank*chunk+j], cols, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
	}
	*/
	//====================================================================
}


