// Matrix addition program MatrixMult.cu, Barry Wilkinson, Dec. 28, 2010.
#include <stdio.h>
#include <cuda.h>
#include <stdlib.h>

__global__ void gpu_matrixmult(int *gpu_a, int *gpu_b, int *gpu_c, int N) {

	int k, sum = 0;
	int col = threadIdx.x + blockDim.x * blockIdx.x; 
	int row = threadIdx.y + blockDim.y * blockIdx.y;

      if (col < N && row < N) {
				for (k = 0; k < N; k++) 
		        		sum += gpu_a[row * N + k] * gpu_b[k * N + col];
				gpu_c[row * N + col] = sum;
			}

}

void cpu_matrixmult(int *cpu_a, int *cpu_b, int *cpu_c, int N) {
	int row, col, k, sum;

	for (row =0; row < N; row++)   				// row of a
		for (col =0; col < N; col++) {				// column of b
			sum = 0;
			for(k = 0; k < N; k++) 
          			sum += cpu_a[row * N + k] * cpu_b[k * N + col];
			cpu_c[row * N + col] = sum;
			//d[row * N + col] = gpu_c[row *N + col];
		}
}


int main(int argc, char *argv[])  {
	int i, j; 							// loop counters
	int Grid_Dim_x=1, Grid_Dim_y=1;		//Grid structure values
	int Block_Dim_x=1, Block_Dim_y=1;		//Block structure values
	int noThreads_x, noThreads_y;			// number of threads available in device, each dimension
	int noThreads_block;					// number of threads in a block
	int N = 10;  						// size of array in each dimension
	int B;
	int T;
	int *a,*b,*c,*d;
	int *dev_a, *dev_b, *dev_c;
	int size;							// number of bytes in arrays
	cudaEvent_t start, stop;     				// using cuda events to measure time
	float elapsed_time_ms;       			// which is applicable for asynchronous code also
	cudaEventCreate(&start);		
	cudaEventCreate(&stop);
int repeat = 1;
while(repeat == 1) {
/* --------------------ENTER INPUT PARAMETERS AND ALLOCATE DATA -----------------------*/
							// keyboard input

	printf("Enter the value for N: ");
	scanf("%d", &N);
//takes in input
	int valid = 0;
	while(valid == 0) {

		printf("Enter the number of blocks: ");
		scanf("%d", &B);

		printf("Enter the number of threads: ");
		scanf("%d", &T);

		if(B > 1024 || T > 1024 || B*T < N*N) {
			printf("Invlaid input entered.");
		} else {
			valid = 1;
			Grid_Dim_x = B;
			Block_Dim_x = T;		//puts the size of blocks and thread in for the dim3
		}
	}
	
	dim3 Grid(Grid_Dim_x, Grid_Dim_x);	//Grid structure
	dim3 Block(Block_Dim_x,Block_Dim_y);	//Block structure, threads/block limited by specific device
	size = N * N * sizeof(int);				// number of bytes in total in arrays

	a = (int*) malloc(size);					//dynamically allocated memory for arrays on host
	b = (int*) malloc(size);
	c = (int*) malloc(size);					// results from GPU
	d = (int*) malloc(size);				// results from CPU
							// load arrays with some numbers

		int row, col;
		srand(2);
		for(row=0; row < N; row++) { // load arrays with some numbers
			for(col=0; col < N; col++) {
				a[row * N + col] = rand() % 10;
				b[row * N + col] = rand() % 10; 
			}
		}

	cudaMalloc((void**)&dev_a, size);			// allocate memory on device
	cudaMalloc((void**)&dev_b, size);
	cudaMalloc((void**)&dev_c, size);

	cudaMemcpy(dev_a, a , size ,cudaMemcpyHostToDevice);
	cudaMemcpy(dev_b, b , size ,cudaMemcpyHostToDevice);

	cudaEventRecord(start, 0); 			// here start time, after memcpy

	gpu_matrixmult<<<Grid,Block>>>(dev_a,dev_b,dev_c,N);
	cudaMemcpy(c, dev_c, size , cudaMemcpyDeviceToHost);

	cudaEventRecord(stop, 0);     			// measuse end time
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsed_time_ms, start, stop );

	printf("Time to calculate results on GPU: %f ms.\n", elapsed_time_ms);
	double gpuTime = elapsed_time_ms; 

/* ------------- COMPUTATION DONE ON HOST CPU ----------------------------*/

	cudaEventRecord(start, 0);			// use same timing*

	cpu_matrixmult(a,b,d,N);				// do calculation on host

	cudaEventRecord(stop, 0);     		// measure end time
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsed_time_ms, start, stop );

	printf("Time to calculate results on CPU: %f ms.\n", elapsed_time_ms);  // exe. time
	double cpuTime = elapsed_time_ms;

/* ------------------- check device creates correct results -----------------*/
/*
	int s;
	for(s=0;s<N*N;s++) {
		printf("%d\t", d[s]);
		if(s%N == 0 && s != 0) {
			printf("\n");
		}
	}
*/
//puts out an error is the two matricies are not the same
	printf("\n");
	int error = 0;
	int k;
	for(k=0; k<N*N; k++) {
		if(d[k] != c[k]) {
			error =1;
			break;
		} 
	}

	if(error ==1 ) 
		printf("There is an error.\n");
	else
		printf("Sequential and parallel produce the same results.\n");

	double speedupFactor;
	speedupFactor = cpuTime/gpuTime;
	printf("Speedup Factor: %lf\n", speedupFactor);

	printf("Would you like to repeat? Enter 1 for yes or 0 for no.\n");
	scanf("%d", &repeat);
}
/* --------------------- repeat program  ----------------------------------------*/
 								//  while loop to repeat calc with different parameters
/* --------------  clean up  ---------------------------------------*/
	free(a); free(b); free(c);
	cudaFree(dev_a);
	cudaFree(dev_b);
	cudaFree(dev_c);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	return 0;
}

