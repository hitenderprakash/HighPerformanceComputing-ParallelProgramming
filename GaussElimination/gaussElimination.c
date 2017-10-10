#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void displayMatrix(double **matrix, int rows, int cols);
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
			mat[i][j]=719+rand()%1971;
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
	
	//Compute inverse
	for(i=0;i<mrow;i++){
		double div=mat[i][i];
		if(div==0){printf("\nError: singular matrix");}
		for(j=0;j<3;j++){
			mat[i][j]=mat[i][j]/div;
			invmat[i][j]=invmat[i][j]/div;
		}
		for(j=0;j<mcol;j++){
			if(i!=j){
				double mul=mat[j][i];
				for(k=0;k<3;k++){
					
					if(mat[j][k]==0 && mat[i][k]==0){}
					else{mat[j][k]=mat[j][k]-(mul*mat[i][k]);}
					
					if(invmat[j][k]==0 && invmat[i][k]==0){}
					else{invmat[j][k]=invmat[j][k]-(mul*invmat[i][k]);}
							
				}
			}
		}
	}
	//inverse calculation section ends
	
	//displayMatrix(mat,mrow,mcol);
	//displayMatrix(invmat,mrow,mcol);

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
