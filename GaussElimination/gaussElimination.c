#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "omp.h"

void displayMatrix(double **matrix, int rows, int cols);
void swaprows(double **matrix,int row1, int row2,int cols);
int main(int argc, char *argv[]){
	int mrow=20000;
	int mcol=20000;
	if(argc!=3){printf("\nUsage: prog <rows> <cols>"); exit(0);}
	mrow=atoi(argv[1]);
	mcol=atoi(argv[2]);
	
	int i,j,k;
	double **mat=0;//initialize pointer with null 
	mat=(double **)malloc(mrow*sizeof(double *));
	
	for(i=0;i<mrow;i++){
		mat[i]=(double *)malloc(mcol*sizeof(double));
	}
	//double count=7;
	//initializing the matrix elements with zero 
	for(i=0;i<mrow;i++){
		for(j=0;j<mcol;j++){
			mat[i][j]=7+rand()%11;
		}
	}
	
	double **invmat=0;//initialize pointer with null 
	invmat=(double **)malloc(mrow*sizeof(double *));
	
	for(i=0;i<mrow;i++){
		invmat[i]=(double *)malloc(mcol*sizeof(double));
	}
	//initializing the matrix elements with zero 
	for(i=0;i<mrow;i++){
		for(j=0;j<mcol;j++){
			if(i==j){invmat[i][j]=1;}
			else{invmat[i][j]=0;}
		}
	}
	//displayMatrix(mat,mrow,mcol);
	//displayMatrix(invmat,mrow,mcol);
	
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
			if(singularFlag==1){printf("\nSingular Matrix ! Exiting");exit(0);}
		}
	}
	
	//Compute inverse
	//Parallel implementation

	for(i=0;i<mrow;i++){
		double div=mat[i][i];
		if(div==0){printf("\nError: singular matrix");exit(0);}
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
	//inverse calculation section ends
	displayMatrix(mat,mrow,mcol);
	displayMatrix(invmat,mrow,mcol);

	printf("\nDone\n");
	return 0;
}

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
void swaprows(double **matrix,int row1, int row2,int cols){
	int i;
	for(i=0;i<cols;i++){
		double temp=matrix[row1][i];
		matrix[row1][i]=matrix[row2][i];
		matrix[row2][i]=temp;
	}
}
