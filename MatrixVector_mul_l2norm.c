/*  Introduction to High Performance Computing 
 *  Assignment 1 
 *  Program for Matrix - vector Multiplication and calculating L2 Norm of resultant Vector
 *  Author: Hitender Prakash
 *  Email: hprakash@iu.edu
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

#define pi 3.14159265
void displayMatrix(double **matrix, int rows, int cols);
double l2NormOfColoumnVector(double **vec, int sz);
double ** matrixMaultiplication(double **leftMatrix, int leftMatrix_rows, int leftMatrix_cols, double **rightMatrix, int rightMatrix_cols);

int main()
{	

	//The number of rows and coloumn in the matrix are 183
	int mrow=183;
	int mcol=183;
	
	/*craeting matrix
	 * This section create the 183x183 matrix 
	 * initiaize each element with 0.0
	 */
	double **mat=0;//initialize pointer with null 
	mat=(double **)malloc(mrow*sizeof(double *));
	
	int i,j;
	for(i=0;i<mrow;i++){
		mat[i]=(double *)malloc(mcol*sizeof(double));
	}
	//initializing the matrix elements with zero 
	for(i=0;i<mrow;i++){
		for(j=0;j<mcol;j++){
			mat[i][j]=0.0;
		}
	}
	
	/*Read data from mtx file and fill the entires in the matrix
	 * The mtx file has sparse matrix entries 
	 * this section parse the file and fill the matrix elements
	 * It assumes that the mtx file is placed in the current directory
	 * 
	 * MTX file has following format:
	 * ROW COLOUMN VALUE
	 * 
	 * Also note that the array indexing starts with [1] in mtx file
	 * while the C Language start array indexing with [0]
	 */
	static const char filename[] = "fs_183_1.mtx";
	FILE *file = fopen ( filename, "r" );
	if ( file != NULL )
	{
       char line [200];
       int line_num=0;
       while ( fgets ( line, sizeof line, file ) != NULL ) 
       {
		  line_num++;
		  
		  //ignore the first two lines of .mtx file, the first line is header line and other represents size info
		  if(line_num<3){ continue;} 
          //parse the matrix info from mtx file
          char *tok;
          tok=strtok(line," ");
          int row,col,tok_num=0;
          double val=0.0;
         
          while(tok!=NULL){
			 tok_num++;		 
			 if(tok_num==1){ row=atoi(tok);} //read row number		 
			 else if(tok_num==2){col=atoi(tok);} //read coloumn number		 
			 else if(tok_num==3){ sscanf(tok, "%lf", &val);} //read value
			 tok=strtok (NULL, " ");
		  }
		  //printf("\nRow: %d  Col: %d  Val: %e",row,col,val);
		  //since in mtx file indexing starts from 1 but in c language array index starts from 0
		  mat[row-1][col-1]=val;
		  //uncomment below line if you want to see that correct value has been written for Matrix[row][col]
		  //printf("\nMatrix[%d][%d] = %e",row,col,mat[row-1][col-1]); 
       }
       fclose ( file );
    }
    else
    {
       printf("\nNOT ABLE TO READ FILE. PLEASE CHECK IF FILE EXISTS");
       printf("\nPLEASE KEEP BOTH, EXECUTABLE AND DEPENDENT MTX FILES IN CURRENT WORKING DIRECTORY");
       printf("\nEXITING NOW...");
       exit(0);
      
    }
	
	//display matrix 
	//please note that this is a sparse matrix and you will see a lot of zero entries
	//uncomment the line below if you want to see the matrix created
	
	//displayMatrix(mat, mrow, mcol); 
	
	/*creating the coloumn vector of dimension 183x1
	 * this section creates the coloumn vector of 183x1 dimension
	 * and each element of the vector is calculated as: Vi = sin ((2* pi* i)/182)
	 */
	int vrow=183;
	int vcol=1;
	double fixVal= (2*pi)/182; 
	double **vec=0; //initialize pointer with null 
	vec=(double **)malloc(vrow*sizeof(double *));
	
	for(i=0;i<vrow;i++){
		vec[i]=(double *)malloc(vcol*sizeof(double));
	}
	
	for(i=0;i<vrow;i++){
		for(j=0;j<vcol;j++){
			vec[i][j]=sin(fixVal*i);
		}
	}
	
	//display vector (uncomment below if you want to see the vector created)
	//displayMatrix(vec, vrow, vcol);

	
	/*multiply matrix-vector
	 * 
	 * 
	 */
	double **resMat=0;
	int resRow=mrow;
	int resCol=vcol;
	resMat=matrixMaultiplication(mat, mrow, mcol, vec, vcol);

	//display result matrix (uncomment if you want to see the resultant vector created)
	//displayMatrix(resMat, resRow, resCol);
	
	//calculate l2Norm of result vector
	double l2norm=l2NormOfColoumnVector(resMat,resRow);
	
	//prints l2norm in scientific( exponential format)
	printf("\nL2 Norm of result vector: %e\n",l2norm);	
	return 0;
}

/* Method for calculating L2 Norm of coloumn vector
 * l2 norm is calculated as: sqrt(x1*x1 + x2*x2 + .... + xn*xn)
 * we can also write it in vector form as: sqrt( transpose(Vector)* vector)
 * Arguents: reference to vector, number of rows in vector
 * Output: L2 Norm of vector
 */
double l2NormOfColoumnVector(double **vec, int sz){
	double sqSum=0.0; //double will cause overflow, taking long double
	int i=0;
	for(i=0;i<sz;i++){
		sqSum=sqSum+vec[i][0]*vec[i][0];
	}
	return sqrt(sqSum);
}

/* Method for displaying the matrix in tabular form
 * Arguments: reference to matrix, number of rows, number of coloumns
 * Output: To console, display of matrix elements in tabular form
 */
void displayMatrix(double **matrix, int rows, int cols){
	if(!matrix){return;}
	int i,j;
	for(i=0;i<rows;i++){
		printf("\n");
		for(j=0;j<cols;j++){
			//using %e to print in scientific notation,
			//use %lf for printing plain double format
			printf("\t%e",matrix[i][j]);
		}
	}
	printf("\n");
}

/* Method for multiplying to Matrices
 * Arguments: left matrix reference, rows in left matrix, coloumns in left matrix, right matrix reference, coloumns in right matrix
 * Assumption: Both the matrices qualify for multiplication i.e coloumns in left matrix are equal to rows in right matrix
 * This method does not validate the above condition and need that valid matrices are provided
 * Input arguments must be sanitized before passing to this method
 * Output: reference to resultant matrix
 */
double ** matrixMaultiplication(double **leftMatrix, int leftMatrix_rows, int leftMatrix_cols, double **rightMatrix, int rightMatrix_cols){
	double **resMat=0;
	int resRow=leftMatrix_rows;
	int resCol=rightMatrix_cols;
	resMat=(double **)malloc(resRow*sizeof(double *));
	
	int i,j,k;
	for(i=0;i<resRow;i++){
		resMat[i]=(double *)malloc(resCol*sizeof(double));
	}
	
	//calculate matrix vector multiplication
	for(i=0;i<resRow;i++){
		for(j=0;j<resCol;j++){
			resMat[i][j]=0.0; //initialize element of result matrix so it does not have garbage value
			for(k=0;k<leftMatrix_cols;k++){
				resMat[i][j]=resMat[i][j]+leftMatrix[i][k]*rightMatrix[k][j];
			}
		}
	}
	return resMat;
}
