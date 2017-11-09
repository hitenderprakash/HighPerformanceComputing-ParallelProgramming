#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define COLS 10
#define ROWS 10
#define TEMP 50.
#define DEBUG 1
#define EPS 1e-6
#define I_FIX 5
#define J_FIX 5

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
			else if(i==0){upper=ghostrow_recv_upper[0][j];}
			if(i!=rows-1){lower=old[i+1][j];}
			else if(i==rows-1){upper=ghostrow_recv_lower[0][j];}			
			newmat[i][j]= (0.25*(upper+lower+right+left));
		}
	}

}

void computeJacobi2(double **oldmat, double **newmat, int rows, int cols, double *ghostrow_recv_upper, double *ghostrow_recv_lower,int rank,int chunk){
	
	int startRow=chunk*rank;	
	int i,j;
	for(i=0;i<chunk;i++){
		for(j=0;j<cols;j++){
			int upper=0;
			int lower=0;
			int right=0;
			int left=0;
			if(j!=0){left=old[i+startRow][j-1];}
			if(j!=cols-1){right=old[i+startRow][j+1];}
			if(i+startRow!=0){upper=old[i+startRow-1][j];}
			else if(i+startRow==0){upper=ghostrow_recv_upper[0][j];}
			if(i+startRow!=rows-1){lower=old[i+startRow+1][j];}
			else if(i+startRow==rows-1){upper=ghostrow_recv_lower[0][j];}			
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

void copy2dTo1d(int *dest, int **source,int size, int rank, int chunk){
		int startRow=rank*chunk;
		int endRow=rank*chunk+chunk-1;
		int i,j;
		for(i=startRow;i<=endRow;i++){
			for(j=0;j<size;j++){			
				dest[i*size+j]=source[i][j];
				//*dest=source[i][j];
				//dest++;
			}
		}
}

void copy1dTo2d(int **dest, int *source,int size, int rank, int chunk){
		int startRow=rank*chunk;
		int endRow=rank*chunk+chunk-1;
		int i,j;
		for(i=startRow;i<=endRow;i++){
			for(j=0;j<size;j++){
				dest[i][j]=source[i*size+j];
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
			max_diff= max(max_diff,newMat[i][j]-oldMat[i+startRow][j];
		}
	}
	return max_diff;
}

int main(int argc, char *argv[]) {
	
	MPI_Status status;
	int rows=ROWS;
	int cols=COLS;
	
	double **a_old;
	double **subMat;
	double **subMatCP;
	
	double maxerr=DBL_MAX;
	double global_err;
	double recverr;
	
	double *ghostrow_send_upper=(double*) calloc(cols,sizeof(double));
	double *ghostrow_send_lower=(double*) calloc(cols,sizeof(double));
	double *ghostrow_recv_upper=(double*) calloc(cols,sizeof(double));
	double *ghostrow_recv_lower=(double*) calloc(cols,sizeof(double));
		
    MPI_Init(&argc,&argv);
	int rank,numProc, root = 0;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&numProc);
    
    int chunk=0; //rows to be taken up by single process
	if(rows%numProc==0){
		chunk=(rows/numProc);
	}
	else{
		chunk=(rows/numProc)+1;
	}
	int startRow=rank*chunk;
	int endRow=min(rank*chunk+chunk-1,cols-1);
		
	//now need to alloc and init the matrix with perfect size (padding may be required)
	a_old=alloc_matrix(chunk*numProc, cols);
	a_old[I_FIX][J_FIX] = TEMP;
	//initially all processes have the same view of the heat sheet the actual averaging starts withing the loop 
	
	//initialize sub matrix with less number of arrays
	subMat=alloc_matrix(chunk, cols);
	subMatCP=alloc_matrix(chunk, cols);
	while(1){
		//send and receive max error to check if we need next iteration
		if(rank!=0){
			MPI_Send(maxerr, 1, MPI_DOUBLE, root, 0, MPI_COMM_WORLD);
		}
		else{
			MPI_Recv(recverr, 1, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD, &status);
			global_err=max(global_err,recverr);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if(rank==0){MPI_Bcast(global_err, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);}	
		if(global_err<EPS){
			break;
		}
		//otherwise proceed 
		
		//rank 0 sends the orig matrix to all procs
		if(rank==0){
			//make sure fix temp is the same at a given point 
			a_old[I_FIX][J_FIX] = TEMP;
			//and also padded elements are zeroed out
			int i,j;	
			for(i=rows;i<chunk*numProc;i++){
				for(j=0;j<cols;j++){
					a_old[i][j]=0.0;
				}
			}
			//send orig matrix
			for(i=1;i<numProc;i++){
				MPI_Send(a_old[rank*chunk], cols*chunk, MPI_DOUBLE, i, 2, MPI_COMM_WORLD);
			}
		}
		else{
			MPI_Recv(a_old[rank*chunk], cols*chunk, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &status);
		}
		
		
		MPI_Barrier(MPI_COMM_WORLD);
		//first task is to send and receive ghost rows
		
		if(rank==0){
			//calculate what to send
			MPI_Send(ghostrow_send_lower, cols, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD);
			MPI_Recv(ghostrow_recv_lower, cols, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &status);
		}
		else if(rank==numProc-1){
			//calculate what to send
			MPI_Send(ghostrow_send_upper, cols, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
			MPI_Recv(ghostrow_recv_upper, cols, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
		}
		else{
			//calculate what to send
			//calculate what to send
			MPI_Send(ghostrow_send_lower, cols, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD);
			MPI_Recv(ghostrow_recv_lower, cols, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &status);
			MPI_Send(ghostrow_send_upper, cols, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
			MPI_Recv(ghostrow_recv_upper, cols, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
		}
		//when finished start computation
		
		computeJacobi2(a_old, subMatCP,rows,cols,ghostrow_recv_upper,ghostrow_recv_lower,rank,chunk);
		
		//compute the max diff here for every thread
		maxerr=findMaxDiff(a_old,subMatCP,cols,rank,chunk);
		
		//copy back the subMatCP to a_old corresponding elements 
		
		copyPartToMatrix(a_old,subMatCP,cols,rank,chunk);
		
		//reconstruct a_old from all procs at proc 0
		if(rank==0){
			//receive from all
			for(i=1;i<numProc;i++){
				MPI_Recv(a_old[rank*chunk], cols*chunk, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &status);
			}
			a_old[I_FIX][J_FIX] = TEMP;
		}
		else{
			//send to 0
			MPI_Send(a_old[rank*chunk], cols*chunk, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		//display
		print_matrix(a_old,rows,cols);
		
	}
	if(rank==0){
		//final display
		printf("\nFinal Matrix: \n");
		print_matrix(a_old,rows,cols);
	}
    
   MPI_Finalize();

    return 0;
}


