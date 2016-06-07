/*******************************************************
File name: primenumber_yahongliu
Author: Yahong Liu
Version: 3.0
Date: 24/may/2016
Description: Sieve of Eratosthenes pthread algorithm
********************************************************/

#include <pthread.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NOT_PRIME 0
#define PRIME 1


// Array for calculation, shared array for every thread
int *array_shared;


// Struct for thread data
typedef struct thread_data
{
	int	thread_id;
	int thread_array_largest;
	int thread_num_threads;
}Thread_data;


// Function declaration
void *thread_function (void *threadarg);


// Main
int main(int argc, char *argv[])
{
	// Largest integer
	int INT_LARGEST = atoi(argv[1]);
	// Number of threads
	int NUM_THREADS = atoi(argv[2]);
	// Check input: INT_LARGEST should be larger than 1 and NUM_THREADS larger than 0.
	if (INT_LARGEST < 2 || NUM_THREADS < 1) {
		printf("Your input is invalid. Please ensure INT_LARGEST is larger than or equal 2 and NUM_THREADS is larger than or equal to 1\n");
		exit(EXIT_FAILURE);
	}
	// Array boundary is INT_LARGEST+1.
	int ARRAY_LARGEST = INT_LARGEST+1;

	// Thread id
	pthread_t threads[NUM_THREADS];
	// Counting variables
	int rc, t, i;
	// Array for storing the thread data structure
	Thread_data thread_data_array[NUM_THREADS];
	// Shared array for storing prime numbers among threads
	array_shared = (int*)malloc( ARRAY_LARGEST * sizeof(int) );
	// set all number to PRIME
	for (int i = 2; i < ARRAY_LARGEST; i++){
		printf("all i: %d\n", i);
		array_shared[i] = PRIME;
	}
	// set even number to NOT_PRIME
	for (i = 2; i < ARRAY_LARGEST+1; i += 2) {
		printf("even i: %d;\n", i);
		array_shared[i] = NOT_PRIME;
	}
	
	// Create thread loop
	for (t=0; t<NUM_THREADS; t++) {
		// set thread_data structure
		thread_data_array[t].thread_id = t;
		thread_data_array[t].thread_array_largest = ARRAY_LARGEST;
		thread_data_array[t].thread_num_threads = NUM_THREADS;

		// create a thread
		printf("Creating thread %d\n", t);
		rc = pthread_create(&(threads[t]), NULL, thread_function, (void *)&thread_data_array[t]);
		if (rc) {
			printf("Thread create failed; return code from pthread_create() is %d\n", rc);
			exit(EXIT_FAILURE);
		}
	}

	printf("Waiting for thread to finish...\n");

	// Join thread loop
	for (t=0; t<NUM_THREADS; t++) {
		// join a thread
		rc = pthread_join(threads[t], NULL);
		if (rc) {
			printf("Thread join failed; return code from pthread_create() is %d\n", rc);
			exit(EXIT_FAILURE);
		}
	}

	printf("Finished searching prime numbers.\n");

	// Special case: number 2 is eliminated as the even number, we need to set it PRIME
	if(array_shared[2] == NOT_PRIME){
		array_shared[2] = PRIME;
	}

	// Open file and output the prime numbers
	FILE *FP;
	char file[] = "primeoutput";
	FP = fopen(file,"w");
	if (!FP){
		fprintf(stderr, "%s cannot be opened!!!\n", file);
	}
	// print prime number to stdout and the file
	for (i = 0; i<ARRAY_LARGEST; i++) {
		if (array_shared[i] == PRIME) {
			printf("%d\n", i);
			fprintf(FP, "%d\n", i);
		}
	}

	// free array_shared and exit program
	free(array_shared);
	exit(EXIT_SUCCESS);
}


// Function implementation
void *thread_function (void *threadarg)
{
	// Thread data variables
	int id;
	int array_largest;
	int num_threads;
	struct thread_data *my_data;

	// Thread boundary and array variables
	int min_boundary, max_boundary, array_size;
	int* array_thread;
	int i, j, k;

	// Parse data from main thread
	my_data = threadarg;
	id = my_data->thread_id;
	array_largest = my_data->thread_array_largest;
	num_threads = my_data->thread_num_threads;

	// For testing: let each thread work in sequence
	// sleep(id+1);

	// Determine each thread boundary
	min_boundary = id * (array_largest/num_threads+1) - ( id<(array_largest % num_threads)?0:1 ) * (id - array_largest % num_threads);
	max_boundary = (id+1) * (array_largest/num_threads+1) - ( id<(array_largest % num_threads)?0:1 ) * (id - array_largest % num_threads + 1) - 1;
	// Eliminate 0 and 1 when min_boundary is smaller than 2
	min_boundary = (min_boundary < 2) ? 2 : min_boundary;
	// Determin array size for each thread
	array_size = max_boundary - min_boundary + 1;
	// Eliminate situation that max_boundary is smaller than 2 or array size is smaller than 1
	if (max_boundary < 2 || array_size < 1){
		pthread_exit(NULL);
	}
	// Determine checking boundary for the Sieve of Eratosthenes
	int checkEnd = sqrt( (double) max_boundary ) + 1;
	printf("thread_id: %d; min_boundary: %d; max_boundary: %d; array_size: %d; checkEnd: %d\n", id, min_boundary, max_boundary, array_size, checkEnd);
	
	// loop for checking divisors, ignore even numbers
	for(i = 3; i < checkEnd; i += 2)
	{
		printf("divisor number loop, i: %d\n", i);
		// if i is composite number, go to next loop
		// if(i >= min_boundary && array_shared[i] == NOT_PRIME){
		if(array_shared[i] == NOT_PRIME){
			printf("divisor is composite number, i: %d\n", i);
			continue;
		}
		
		// check every 2 number from the first odd number in this thread, also ignore even numbers
		for (j = ((min_boundary%2==0)?(min_boundary+1):min_boundary); j < max_boundary+1; j += 2)
		{
			printf("array number loop, j: %d\n", j);
			// check if the array number is composite number
			if (j!=i && j%i == 0)
			{
				printf("j%%i==0, i: %d; j: %d\n", i, j);
				array_shared[j] = NOT_PRIME;
			}
		}
	}

	// exit thread
	pthread_exit(NULL);
}
