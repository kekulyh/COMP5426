#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/time.h>
#include <mpi.h>

// define BOOL type
#ifndef BOOL
#define BOOL int
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

/*******************************************************
File name: redblue_yahongliu
Author: Yahong Liu
Version: 4.0
Date: 15/April/2016
Description: red blue computation parallel algorithm
********************************************************/

/*
 white = 0, red = 1, blue = 2,
 red or blue just moved in = 3 and
 red or blue (in the first row or column) just moved out = 4
 */

// Function Declaration
void create_grid(int ***grid, int gridWidth, int gridHeight);

void board_init(int ***grid, int gridWidth, int gridHeight);

void display_grid(int **grid, int gridWidth, int gridHeight);

void move_grid_sequential(int ***grid, int n, int *n_itrs, int MAX_ITRS, BOOL finished, int t, int c, int redcount, int bluecount);

int check_tile(int **grid, int redcount, int bluecount, int n, int t, int c, BOOL finished, int n_itrs );

// Functions Implementation
// create dynamic two-dimensional array
void create_grid(int ***grid, int gridWidth, int gridHeight) {
	
	int i;
	
	*grid = malloc( sizeof( int* ) * gridHeight ); // create a column
	
	for (i = 0; i < gridHeight; i++){
		(*grid)[i]=(int *)malloc(gridWidth*sizeof(int)); // create rows for previous column
	}
}

// initialize the grid or block with colored squares
void board_init(int ***grid, int gridWidth, int gridHeight) {
	
	// linux struct to get microsecond
	struct timeval tv;
	// get time
	gettimeofday(&tv, NULL);
	// set srand function for randomly generategrid
	// tv_usec: microseconds; tv_sec: seconds.
	srand(tv.tv_usec);

	int i, j, k;
	
	// Assign colors for each square with 1/3 probability
	for(i = 0; i < gridHeight; i++)
	{
		for(j = 0; j < gridWidth; j++){
			
			k = rand() % 100;
			if ( k < 33 ) {
				(*grid)[i][j] = 0;
			} else if ( k < 66 ){
				(*grid)[i][j] = 1;
			} else {
				(*grid)[i][j] = 2;
			}
		}
	}
}

// display the grid, must be annotated when grid size is too large
void display_grid(int **grid, int gridWidth, int gridHeight) {
	
	int i, j;
	
	for(i = 0; i < gridHeight; i++)
	{
		for(j = 0; j < gridWidth; j++)
		{
			printf("%d ", grid[i][j]);
		}
		printf("\n");
	}
    
}

// Sequential function for moving color squares
void move_grid_sequential(int ***grid, int n, int *n_itrs, int MAX_ITRS, BOOL finished, int t, int c, int redcount, int bluecount) {
	
	int i, j;
	
	while ( !finished && (*n_itrs) < MAX_ITRS ){
		(*n_itrs)++;
		
		/* red color movement */
		for (i = 0; i < n; i++){
			if ((*grid)[i][0] == 1 && (*grid)[i][1] == 0){
				(*grid)[i][0] = 4;
				(*grid)[i][1] = 3;
			}
			for (j = 1; j < n; j++){
				if ((*grid)[i][j] == 1 && ((*grid)[i][(j+1)%n] == 0)){
					(*grid)[i][j] = 0;
					(*grid)[i][(j+1)%n] = 3;
				} else if ((*grid)[i][j] == 3){
					(*grid)[i][j] = 1;
				}
			}
			if ((*grid)[i][0] == 3){
				(*grid)[i][0] = 1;
			} else if ((*grid)[i][0] == 4){
				(*grid)[i][0] = 0;
			}
			// printf("Red move:\n");
			// display_grid((*grid), n);
		}
		
		/* blue color movement */
		for (j = 0; j < n; j++){
			if ((*grid)[0][j] == 2 && (*grid)[1][j] == 0){
				(*grid)[0][j] = 4;
				(*grid)[1][j] = 3;
			}
			for (i = 1; i < n; i++){
				if ((*grid)[i][j] == 2 && (*grid)[(i+1)%n][j]==0){
					(*grid)[i][j] = 0;
					(*grid)[(i+1)%n][j] = 3;
				}
				else if ((*grid)[i][j] == 3){
					(*grid)[i][j] = 2;
				}
			}
			if ((*grid)[0][j] == 3){
				(*grid)[0][j] = 2;
			}
			else if ((*grid)[0][j] == 4){
				(*grid)[0][j] = 0;
			}
			// printf("Blue move:\n");
			// display_grid((*grid), n);
		}
		
		// invoke check_tile to get 'finished' value(true or false)
		finished = check_tile((*grid), redcount, bluecount, n, t, c, finished, (*n_itrs));
		
	}
	
	// maximum iteration has reached, remind user to input a larger iteration number
	if ( finished == FALSE )
	{
		printf("Iteration number %d reached, pleaase input a larger iteration number.\n", (*n_itrs));
	}
	
}

// function for checking colored squares in each tile
int check_tile(int **grid, int redcount, int bluecount, int n, int t, int c, BOOL finished, int n_itrs ) {
	
	int i, j, k, l;

	int tile_size = n/t;
	
	for (i = 0; i < t; i++){
		for (j = 0; j < t; j++){ //loop in tiles in each row and colomn
			
			redcount = 0;
			bluecount = 0;
			
			for (k = i * tile_size; k < tile_size * ( i + 1 ); k++) {
				for (l = j * tile_size; l < tile_size * ( j + 1 ); l++) { //loop in each tile
					
					// count the colored squares
					switch (grid[k][l]){
						case 0:
							break;
							
						case 1:
							redcount += 1;
							break;
							
						case 2:
							bluecount += 1;
							break;
					}
					
				}
			}
			// print the counts and percentage of this tile
			printf("-------------------------------------\n");
			printf("[%d, %d] tile:\n", i, j);
			printf("redcount: %d\n", redcount);
			printf("bluecount: %d\n", bluecount);
			printf("red color percentage: %f %%\n", ( (redcount*100.00) / (tile_size*tile_size) ) );
			printf("blue color percentage: %f %%\n", ( (bluecount*100.00) / (tile_size*tile_size)) );
			
			// if counts of one color has reached the threshold, then terminate the loop
			if ( ( ( (redcount*100.00) / (tile_size*tile_size) ) > c ) || ( ( (bluecount*100.00) / (tile_size*tile_size) ) > c ) )
			{
				printf("-------------------------------------\n");
				printf("Stop computing, [%d, %d] tile reaches thresh %d %%, iteration: %d.\n", i, j, c, n_itrs);
				finished = TRUE;
				break;
			}
			
		}
		
		if (finished)
		{
			break;
		}
		
	}
	
	// return 'finished' value to its invoker
	return finished;
}


int main(int argc, char *argv[]) {
	
	// MPI variables
	// myid and numprocs are for group use
	// myid_world and numprocs_world are for global use
	int myid, numprocs, myid_world, numprocs_world;

	int source, destination;
	
	MPI_Status status;
	MPI_Init(&argc,&argv); /* Used to send the command line argumenys to all procs */
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs_world); /* Initializes the number of procs in the group specified by mpirun */
	MPI_Comm_rank(MPI_COMM_WORLD,&myid_world);	 /* Initialize the rank of this process in the group */
	
	// Define variables for input error check
	char* errorSize, *errorTile, *errorThreshold, *errorIters; //pointer to character with wrong input
	errno = 0; //number of last error
	
	// Get input arguments
	int n = strtol(argv[1], &errorSize, 10); //grid size
	int t = strtol(argv[2], &errorTile, 10); //tile size: block numbers of each row or column
	int c = strtol(argv[3], &errorThreshold, 10); //terminating threshold
	int MAX_ITRS = strtol(argv[4], &errorIters, 10); //maximum number of iterations
	
	// Determine whether input arguments are correct
	if ((errno == ERANGE) || (errno != 0) || ((*errorSize != 0) || (*errorTile != 0) || (*errorThreshold != 0) || (*errorIters != 0))) {
		perror("Error input arguments");
		exit(EXIT_FAILURE);
	}
	
	// Define variables for grid computation
	int **grid; //two-dimensional array
	BOOL finished = FALSE; // 'finished' flag
	int n_itrs = 0; //current iteration
	int redcount, bluecount;
	int tile_size = n/t; // how many squares in each row or column of each tile
	

	printf("MPI start...\n");
	
	// MPI start
	// Only one processor, sequential program
	if (numprocs_world==1)
	{
		
		printf("Sequential Program\n");
		printf("Numprocs: %i; Process: %i.\n", numprocs_world, myid_world);

		printf("-------------------------------------\n");
		
		// Start Computing
		// Create grid
		printf("Creating grid...\n");
		create_grid(&grid, n, n); // pass '&p' to function，so inside the function, use '*(&p)=p' to modify the global arrray
		printf("Creating sucess\n");
		
		printf("-------------------------------------\n");
		
		// Initialize grid
		printf("Initializing grid...\n");
		board_init(&grid, n, n);
		printf("Initializing sucess\n");
		
		printf("-------------------------------------\n");
		
		// Displaying grid
		// printf("Displaying original grid:\n");
		// display_grid(grid, n, n); // pass 'p' to function，so inside the function, use 'p' to access the arrray. But the function cannot modify the global array
		
		printf("-------------------------------------\n");

		// Move grid
		printf("Moving grid:\n");
		move_grid_sequential(&grid ,n, &n_itrs, MAX_ITRS, finished, t, c, redcount, bluecount);
		
		printf("-------------------------------------\n");
		
		// free the dynamic array
		free(grid);
		
		printf("Program finished\n");
	}
	else // multiple processors, parallel program
	{
		printf("Parallel Program\n");

		int i, j, k, l;

		// Define group id of processors
		int group_id = myid_world / t;

		// Define new MPI_Comm
		MPI_Comm MPI_COMM_GROUP;
		/* split the original communicator based on the color and use the original rank for ordering */
		MPI_Comm_split(MPI_COMM_WORLD, group_id, myid, &MPI_COMM_GROUP);

		/* get the size and rank in the new communicator */
		MPI_Comm_size(MPI_COMM_GROUP, &numprocs);
		MPI_Comm_rank(MPI_COMM_GROUP, &myid);

		/* print out original id, group id in the group */
		printf("my original id is %i, group id is %i in group %i\n", myid_world, myid, group_id);
		printf("-------------------------------------\n");

		// Define the id of previous and next processor
		destination = (myid + 1) % numprocs; //next processor
		source = (myid -1 + numprocs) % numprocs; //previousprocessor
		
		// if not working processors, block them
		if(group_id != 0){

			MPI_Barrier(MPI_COMM_GROUP);

		}
		else{ // if working processors, start computation

		
			// Define edges and size of this processor
			// 我们先假设每个processor都有多分配的一个tile行，这样计算下来，会超出grid的边界
			// myid * (t / numprocs + 1) * tile_size  是假设的该processor第一行
			// (myid + 1)* (t / numprocs + 1) * tile_size  是假设的该processor最后一行
			// 然而只有满足 myid<(t % numprocs) 条件，该processor才会实际上有多分配的一个tile行
			// 所以，对于不满足该条件的processor，减去相应的多分配的，则可以得到正确的边界

			// int PROC_MIN = myid * (t/numprocs+1) * tile_size  - ( myid<(t % numprocs)?0:1 ) * (myid - t % numprocs) * tile_size;
			// int PROC_MAX = (myid+1) * (t/numprocs+1) * tile_size - ( myid<(t % numprocs)?0:1 ) * (myid - t % numprocs + 1) * tile_size - 1 ;
			int PROC_SIZE = ( t/numprocs - ( myid<(t % numprocs)?0:1 ) + 1 ) * tile_size;

			// print processor edges for testing
			printf("Numprocs: %i; Processor: %i; PROC_SIZE: %d\n", numprocs, myid, PROC_SIZE);
			printf("-------------------------------------\n");

			// Create grid
			// Each processor create its own blocks
			create_grid(&grid, n, PROC_SIZE);

			// Initialize grid
			board_init(&grid, n, PROC_SIZE);
			
			// Display grid
			// display_grid(grid, n, PROC_SIZE);

			// Two buffers for each processor
			int *bufferUp = malloc( sizeof(int) * n );
			int *bufferDown = malloc( sizeof(int) * n );

			// Computation starts
			while ( !finished && n_itrs < MAX_ITRS ){
				n_itrs++;
				
				// 1. move red squares in each processor
				/* red color movement */
				for (i =  0; i < PROC_SIZE; i++){
					if (grid[i][0] == 1 && grid[i][1] == 0){
						grid[i][0] = 4;
						grid[i][1] = 3;
					}
					for (j = 1; j < n; j++){
						if (grid[i][j] == 1 && (grid[i][(j+1)%n] == 0)){
							grid[i][j] = 0;
							grid[i][(j+1)%n] = 3;
						} else if (grid[i][j] == 3){
							grid[i][j] = 1;
						}
					}
					if (grid[i][0] == 3){
						grid[i][0] = 1;
					} else if (grid[i][0] == 4){
						grid[i][0] = 0;
					}
					// printf("Red move:\n");
					// display_grid(grid, n);
				}

				// 2. send maximum edge to bufferUp of next processor, and get maximum edge of previous processor into its own bufferUp
				MPI_Sendrecv(grid[PROC_SIZE-1], n, MPI_INT, destination, 1, bufferUp, n, MPI_INT, source, 1, MPI_COMM_GROUP, &status);

				// 3. send minimum edge to bufferDown of previous processor, and get minimum edge of next proccessor into its own bufferDown
				MPI_Sendrecv(grid[0], n, MPI_INT, source, 1, bufferDown, n, MPI_INT, destination, 1, MPI_COMM_GROUP, &status);

				// 4. compare the minimum edge with bufferUp and move blue squares
				for (k = 0; k < n; k++)
				{
					if (bufferUp[k] == 2 && grid[0][k]==0){
						bufferUp[k] = 4;
						grid[0][k] = 3;
					}
				}

				// 5. move the blue squares of processor
				/* blue color movement */
				for (j = 0; j < n; j++){
					if (grid[0][j] == 2 && grid[1][j] == 0){
						grid[0][j] = 4;
						grid[1][j] = 3;
					}
					// The maximum edge cannot move to its next rows otherwise array will get out of range. So the maximum cycle is 'PROC_SIZE - 1'
					for (i = 0 ; i < PROC_SIZE - 1; i++){
						if (grid[i][j] == 2 && grid[(i+1)%n][j]==0){
							grid[i][j] = 0;
							// 这行的(i+1)%n应该为(i+1)%PROC_SIZE
							grid[(i+1)%n][j] = 3;
						}
						else if (grid[i][j] == 3){
							grid[i][j] = 2;
						}
					}
					if (grid[0][j] == 3){
						grid[0][j] = 2;
					}
					else if (grid[0][j] == 4){
						grid[0][j] = 0;
					}
					// printf("Blue move:\n");
					// display_grid(grid, n);
				}

				// 6. compare the maximum edge with bufferDown and move blue squares
				for (l = 0; l < n; l++)
				{
					if (bufferDown[l] == 0 && grid[PROC_SIZE-1][l]==2){
						bufferDown[l] = 3;
						grid[PROC_SIZE-1][l] = 4;
					}
				}

				// 7. check colored squares in each tile of processor
				int p, q, r, s;
				// loop in tiles in each row and colomn
				for (p = 0; p < ( PROC_SIZE/tile_size ); p++){ // each column has p blocks
					for (q = 0; q < t; q++){ // each row has q blocks

						redcount = 0;
						bluecount = 0;

						for (r = p * tile_size; r < tile_size * ( p + 1 ); r++) {
							for (s = q * tile_size; s < tile_size * ( q + 1 ); s++) { //loop in each tile

								switch (grid[r][s]){
									case 0:
										break;
									case 1:
										redcount += 1;
										break;
									case 2:
										bluecount += 1;
										break;
								}

							}
						}

						// print the counts and percentage of this tile
						printf("-------------------------------------\n");
						printf("processor: %d, [%d, %d] tile:\n", myid, p, q);
						printf("processor: %d, redcount: %d\n", myid, redcount);
						printf("processor: %d, bluecount: %d\n", myid, bluecount);
						printf("processor: %d, red color percentage: %f %%\n", myid, ( (redcount*100.00) / (tile_size*tile_size) ) );
						printf("processor: %d, blue color percentage: %f %%\n", myid, ( (bluecount*100.00) / (tile_size*tile_size)) );

						// if counts of one color has reached the threshold, then terminate the loop
						if ( ( ( (redcount*100.00) / (tile_size*tile_size) ) > c ) || ( ( (bluecount*100.00) / (tile_size*tile_size) ) > c ) )
						{
							printf("-------------------------------------\n");
							printf("Processor %d Stop computing, [%d, %d] tile reaches thresh %d %%, iteration: %d.\n", myid, p, q, c, n_itrs);
							finished = TRUE;
							break;
						}

					}

					if (finished)
					{
						break;
					}

				}
				
				// 10. Reduce all the 'finished' flag and do the logical OR operation.
				// Thus, if any processor has set 'finished' flag to TRUE, the result of logical OR will be TRUE.
				// Then, broadcast this 'finished' flag to each processors to terminate programs.
				MPI_Allreduce(&finished, &finished, 1, MPI_INT, MPI_LOR, MPI_COMM_GROUP);
				
			} // while loop end

			// maximum iteration has reached, remind user to input a larger iteration number
			if ( finished == FALSE )
			{
				printf("Processor %d Iteration number %d reached, pleaase choose a larger iteration number.\n", myid, n_itrs);
			}

			printf("-------------------------------------\n");

			// free the dynamic array of each processor
			free(grid);
			free(bufferUp);
			free(bufferDown);

			printf("Processor %d Program finished\n", myid);


		} //working processor end


	} // parallel program end

	// MPI finalize
	MPI_Finalize();

	return 0;
}