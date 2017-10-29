/*
 * Program to sort integers - Parallel Quicksort using MPI
 * Hitender Prakash
 * Approach: process with rank 0 initialize an array with rand()
 * It divide the array in chunks and scatter within all processes
 * each process sort its chunk using quick sort (sequential)
 * Process with rank 0 gathers the chunks back to original array
 * which is sorted in chunks
 * Process with rank 0, now sequentially merge the sorted chunks to a new array
 * 
 * version 2: cosmetic changes and removing unwanted code
 * Date: Oct 28, 2017
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define min(a,b) (((a)<(b))?(a):(b))
#define N 21

//swap two elements in array
void swap(int *arr, int i, int j){
	if(!arr || i==j){return;}
	int temp=arr[i];
	arr[i]=arr[j];
	arr[j]=temp;
}

//display an array between two indices
void displayArray(int *arr, int l, int h){
	if(!arr){return;}
	int i;
	printf("\n");
	for(i=l;i<=h;i++){ printf("%d ",arr[i]);}
	printf("\n");
}
//=====================================================================
// classical quick sort implementation with partition() and quickSort()
int partition(int *arr, int l, int h){
	if(l>h){return -1;}
	if(l==h){return l;}
	//pivot is the last element of the range i.e at 'h' index
	int i=l-1;
	int j=l;
	while(j<h){
		if(arr[j]<= arr[h]){swap(arr, ++i,j);}
		j++;
	}
	swap(arr,++i,h);
	return i;
}

void quickSort(int *arr, int l, int h){
	if(arr && l<h){
		int p=partition(arr, l, h);
		quickSort(arr,l,p-1);
		quickSort(arr,p+1,h);
	}
}
//======================================================================

//Method to merge the p sorted chunks into a single array 
void SequentialMerge(int *oldarr,int *newarr,int chunk, int numProc){
	int itr=0;
	int val;
	//array which keeps track of smallest element in every chunk
	int *procRanges=(int*)malloc(numProc*sizeof(int));
	int i;
	for(i=0;i<numProc;i++){
		procRanges[i]=i*chunk;
	}
	while(itr<(chunk*numProc)){
		val=9999; //initialize val with a large number like INT_MAX
		int incProc;
		for(i=0;i<numProc;i++){
			if(procRanges[i] < (i+1)*chunk){
				if(val >= oldarr[procRanges[i]]){
					val=oldarr[procRanges[i]];
					incProc=i;
				}
			}
		}
		//increment the counter corresponding to chunk from where smallest element was taken
		procRanges[incProc]=procRanges[incProc]+1;
		//write smallest element to new sorted global array
		newarr[itr]=val;
		itr++;
	}	
}

int main(int argc, char **argv)
{

	int size=N; //number of elements to be created and sorted
	int i=0;
	int *origarr;//contains original array
	int *newarr; //sorted return array
	int *sct; //scattered chunk
	
	
	MPI_Init(&argc,&argv);
	int rank,numProc, root = 0;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&numProc);
	
	int chunk=0; //chunk size to be taken up by single process
	if(size%numProc==0){
		chunk=(size/numProc);
	}
	else{
		chunk=(size/numProc)+1;
	}
	sct = (int *)malloc(chunk*sizeof(int));
	if(rank==0){
		newarr=(int *)malloc((chunk*numProc)*sizeof(int));
		origarr=(int *)malloc((chunk*numProc)*sizeof(int));
		for(i=0;i<size;i++){
			origarr[i]=(int)rand()%20;
		}
		//for scattering we want each process to have equal elements in chunk 
		//if array can not be evenly divided we pad it with extra elements
		//padding elements should be INT_MAX so that they remain in the end of the sorted array 
		//and we can ignore them
		for(i=size;i<(chunk*numProc);i++){
			origarr[i]=9999;
		}
		printf("\n");
		printf("\nThe elements to be sorted are:");
		printf("\n=====================================================");
		//displayArray(origarr, 0, (chunk*numProc)-1);
		displayArray(origarr, 0, size-1);	//ignore padded elements			
	}
	
	MPI_Scatter(origarr,chunk,MPI_INT,sct,chunk,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	//sort each chunk by quick sort sequentially
	
	quickSort(sct,0,chunk-1);
	//displayArray(sct, 0, chunk-1);
	
	MPI_Gather(sct,chunk, MPI_INT,origarr, chunk, MPI_INT,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(rank==0){
		//print original array after sorted in chunks
		//displayArray(origarr, 0, (chunk*numProc)-1);
		
		//merge sequentially sorted chunks in a sorted array
		SequentialMerge(origarr,newarr,chunk, numProc);
		
		//print final sorted array
		//displayArray(newarr, 0, (chunk*numProc)-1);
		printf("\nThe elements after sorting are:");
		printf("\n=====================================================");
		displayArray(newarr, 0, size-1); //ignore padded element
		printf("\n\n");
	}

	MPI_Finalize();
	//free(origarr);
	//free(newarr);
	return 0;
}
