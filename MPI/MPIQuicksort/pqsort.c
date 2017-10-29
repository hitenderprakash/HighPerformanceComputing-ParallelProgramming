/*
 * Hitender Prakash
 * version 1
 * Program to sort an array in parallel using quick sort
 * 
 */ 
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define min(a,b) (((a)<(b))?(a):(b))
#define N 21
typedef struct {
	int init;
	int end;
} range;
void swap(int *arr, int i, int j){
	if(!arr || i==j){return;}
	int temp=arr[i];
	arr[i]=arr[j];
	arr[j]=temp;
}
void displayArray(int *arr, int l, int h){
	if(!arr){return;}
	int i;
	printf("\n");
	for(i=l;i<=h;i++){ printf("%d ",arr[i]);}
	printf("\n");
}

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
/*
int *SequentialMerge(int *oldarr,int size, int numProc, range *procRange){
	int *newarr=(int *)malloc(size*sizeof(int));
	int itr=0;
	while(itr<size){
		int val=oldarr[0];
		int incProc;
		int i;
		for(i=0;i<numProc;i++){
			if(procRange[i].init <= procRange[i].end){
				if(val >= oldarr[procRange[i].init]){
					val=oldarr[procRange[i].init];
					incProc=i;
				}
			}
		}
		procRange[incProc].init=procRange[incProc].init+1;
		newarr[itr]=val;
		itr++;
	}	
	return newarr;
}*/
void SequentialMerge(int *oldarr,int *newarr,int chunk, int numProc){
	int itr=0;
	int val;
	int *procRanges=(int*)malloc(numProc*sizeof(int));
	int i;
	for(i=0;i<numProc;i++){
		procRanges[i]=i*chunk;
	}
	while(itr<(chunk*numProc)){
		val=9999;
		int incProc;
		for(i=0;i<numProc;i++){
			if(procRanges[i] < (i+1)*chunk){
				if(val >= oldarr[procRanges[i]]){
					val=oldarr[procRanges[i]];
					incProc=i;
				}
			}
		}
		procRanges[incProc]=procRanges[incProc]+1;
		newarr[itr]=val;

		itr++;
	}	
}

int main(int argc, char **argv)
{
	//create a random array 
	int size=N;
	int i=0;
	int *origarr;
	int *newarr;
	int *sct;
	
	
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
			//newarr[i]=1;
		}
		for(i=size;i<(chunk*numProc);i++){
			origarr[i]=9999;
			//newarr[i]=1;
		}
		displayArray(origarr, 0, (chunk*numProc)-1);	
		//displayArray(newarr, 0, (chunk*numProc)-1);
			
	}
	
	MPI_Scatter(origarr,chunk,MPI_INT,sct,chunk,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	//printf("\nMyRank: %d",rank);
	quickSort(sct,0,chunk-1);
	//displayArray(sct, 0, chunk-1);
	MPI_Gather(sct,chunk, MPI_INT,origarr, chunk, MPI_INT,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0){
		//newarr=(int *)malloc((chunk*numProc)*sizeof(int));
		SequentialMerge(origarr,newarr,chunk, numProc);
		displayArray(newarr, 0, (chunk*numProc)-1);
	}
	/*
	//section to divide the task between processes
	int chunk=0; //chunk size to be taken up by single process
	if(size%numProc==0){
		chunk=(size/numProc);
	}
	else{
		chunk=(size/numProc)+1;
	}
  
	//debug
	//printf("\nChunk Size: %d",chunk);
  
  
	//decide loop boundary to be executed by process 
	int startLoop=chunk*rank;
	int endLoop=min(startLoop+chunk-1,size-1);
	//printf("\nI am rank: %d  I am picking from ind: %d to ind: %d",rank, startLoop, endLoop);
	quickSort(origarr, startLoop, endLoop);
	//printf("\nBang!!");
	procRange[rank].init=startLoop;
	procRange[rank].end=endLoop;
	//debug info
	printf("\nAfter Local Sort:");
	displayArray(origarr, startLoop, endLoop);
	
	MPI_Barrier(MPI_COMM_WORLD);
	printf("\nProc range Table info");
	for(i=0;i<numProc;i++){
		printf("\n %d %d", procRange[rank].init,procRange[rank].end);
	}
	//displayArray(origarr, 0, size-1);
	//int *newarr = SequentialMerge(origarr,size, numProc, procRange);
	//displayArray(newarr, 0, size-1);
	*/
	MPI_Finalize();
	//free(origarr);
	//free(newarr);
	return 0;
}
