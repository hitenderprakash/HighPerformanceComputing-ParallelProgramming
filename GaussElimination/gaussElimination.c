/*
 * Parallel (OpenMP) program for computing Matrix Inverse
 * Algorithm: Guass Elimination (row reduction)
 * Author: Hitender Prakash
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "omp.h"

//function prototypes

//write a string to a file
void write_to_file(char *message,char *filename, char *file_mode);

//write a matrix to the file
void write_matrix_to_file(double** matrix,int rows, int cols, char *filename, char *file_mode);

//display matrix to console
void displayMatrix(double **matrix, int rows, int cols);

//swap two rows in a matrix
void swaprows(double **matrix,int row1, int row2,int cols);

//=======================main programs starts===========================
int main(int argc, char *argv[]){
	int mrow=20000;
	int mcol=20000;
	if(argc!=3){
		printf("\nUsage: prog <rows> <cols>\n"); 
		exit(0);
	}
	//read rows and columns from the arguments 
	mrow=atoi(argv[1]);
	mcol=atoi(argv[2]);
	
	//check if rows and columns are not zero
	if(!mrow || !mcol){
		printf("\nError: Rows or Coloumns cannot be 0 ! Exiting...\n");
		exit(0);
	}
	
	//check if matrix is a square matrix
	
	if(mrow!=mcol){
		printf("\nError: Matrix is not a square matrix! Exiting...\n");
		exit(0);
	}
	
	int i,j,k;
	double **mat=0;//initialize pointer with null 
	mat=(double **)malloc(mrow*sizeof(double *));
	if(!mat){
		printf("\nError: Memory allocation failure ! Exiting...\n");
		exit(0);
	}
	for(i=0;i<mrow;i++){
		mat[i]=(double *)malloc(mcol*sizeof(double));
		if(!mat[i]){
			printf("\nError: Memory allocation failure ! Exiting...\n");
			exit(0);
		}
	}
	//double count=7;
	//initializing the matrix elements with randomly generated values  
	for(i=0;i<mrow;i++){
		for(j=0;j<mcol;j++){
			mat[i][j]=7+rand()%11;
		}
	}
	
	//matrix to store inverse
	double **invmat=0;//initialize pointer with null 
	invmat=(double **)malloc(mrow*sizeof(double *));
	if(!invmat){
		printf("\nError: Memory allocation failure ! Exiting...\n");
		exit(0);
	}
	
	for(i=0;i<mrow;i++){
		invmat[i]=(double *)malloc(mcol*sizeof(double));
		if(!invmat[i]){
			printf("\nError: Memory allocation failure ! Exiting...\n");
			exit(0);
		}
	}
	//initializing the inverse matrix with identity matrix
	for(i=0;i<mrow;i++){
		for(j=0;j<mcol;j++){
			if(i==j){invmat[i][j]=1;}
			else{invmat[i][j]=0;}
		}
	}
	
	//write original matrix data to file 
	write_to_file("Original Input matrix is:\n","results.txt", "w+");
	write_matrix_to_file(mat,mrow,mcol,"results.txt", "a");
	
	//transform if any diagonal element is zero
	int singularFlag=0;
	for(i=0;i<mrow;i++){
		if(mat[i][i]!=0){continue;}
		else{
			singularFlag=1;
			for(j=0;j<mrow;j++){
				if(mat[j][i]!=0 && mat[i][j]!=0){
					singularFlag=0;
					swaprows(mat,i,j,mcol);
					swaprows(invmat,i,j,mcol);
					break;
				}
			}
			if(singularFlag==1){
				write_to_file("\nError: This is a Singular Matrix ! Inverse does not exist!\n","results.txt", "a");
				printf("\nError: This is a Singular Matrix ! Exiting...\n");
				exit(0);
			}
		}
	}
	
	//Compute inverse
	//Parallel implementation

	for(i=0;i<mrow;i++){
		double div=mat[i][i];
		if(div==0){
			write_to_file("\nError: This is a Singular Matrix ! Inverse does not exist!\n","results.txt", "a");
			printf("\nError: This is a Singular Matrix ! Exiting...\n");
			exit(0);
		}
		
		#pragma omp parallel for
		for(j=0;j<mcol;j++){
			mat[i][j]=mat[i][j]/div;
			invmat[i][j]=invmat[i][j]/div;
		}
		#pragma omp barrier
		#pragma omp parallel for shared (i) private(j) private (k)
		for(j=0;j<mrow;j++){
			if(i!=j){
				double mul=mat[j][i];
				for(k=0;k<mcol;k++){
					mat[j][k]=mat[j][k]-(mul*mat[i][k]);
				    invmat[j][k]=invmat[j][k]-(mul*invmat[i][k]);						
				}
			}
		}
		#pragma omp barrier
	}
	
	//write inverse matrix data to file
	write_to_file("\nInverse of the matrix is:\n","results.txt", "a");
	write_matrix_to_file(invmat,mrow,mcol,"results.txt", "a");
	
	//free memory
	for(i=0;i<mrow;i++){
		free(mat[i]);
		free(invmat[i]);
	}
	free(mat);
	free(invmat);
	//===============

	return 0;
}
//main ends

//============================= define functions =======================
/* display matrix to console
 * input:
 * matrix: matrix to be written
 * rows: number of rows in matrix
 * cols: number of columns in matrix
 * Output: void
 */
void displayMatrix(double **matrix, int rows, int cols){
	if(!matrix){return;}
	int i,j;
	for(i=0;i<rows;i++){
		printf("\n");
		for(j=0;j<cols;j++){

			printf("\t%f",matrix[i][j]);
		}
	}
	printf("\n");
}
/* swap two rows in a matrix
 * input:
 * matrix: matrix to be written
 * row1 and row2: rows to be swaped
 * cols: number of columns
 * Output: void
 */
void swaprows(double **matrix,int row1, int row2,int cols){
	int i;
	for(i=0;i<cols;i++){
		double temp=matrix[row1][i];
		matrix[row1][i]=matrix[row2][i];
		matrix[row2][i]=temp;
	}
}
/* write a matrix data to file
 * input:
 * matrix: matrix to be written
 * rows: number of rows in matrix
 * cols: number of columns
 * filename: file to write
 * file_mode: mode to open the file like create or append
 * Output: void
 */
void write_matrix_to_file(double** matrix,int rows, int cols,char *filename, char *file_mode){
	FILE *fout=fopen(filename,file_mode);
	int i,j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++){
            fprintf(fout,"%lf\t", matrix[i][j]);
        }
        fprintf(fout,"\n");
    }
    fclose(fout);
}
/* write a string to a file
 * input:
 * message: string to write
 * filename: file to write
 * file_mode: mode to open the file like create or append
 * Output: void
 */
void write_to_file(char *message,char *filename, char *file_mode){
	FILE *fout=fopen(filename,file_mode);
	fprintf(fout,"%s",message);
	fclose(fout);
}
//============================= functions definition ends =======================
