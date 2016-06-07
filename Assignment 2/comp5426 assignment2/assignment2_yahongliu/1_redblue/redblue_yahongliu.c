/*******************************************************
File name: redblue_yahongliu
Author: Yahong Liu
Version: 4.0
Date: 24/May/2016
Description: red blue computation parallel algorithm
********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/time.h>
#include <mpi.h>
#include <math.h>

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
	int myid, numprocs, numprocs_group, myid_group, myid_world, numprocs_world;

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
	if (numprocs_world == 1)
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
		// printf("-------------------------------------\n");
		
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

		// processor rank larger than t*t is unneeded
		// Define group id of processors
		int group_id = myid_world / ( t * t );

		// Define new MPI_Comm
		// useful group
		MPI_Comm MPI_COMM_GROUP;
		// obseleted group
		MPI_Comm MPI_COMM_GROUP_OBSOLETE;
		
		/* split the original communicator based on the color and use the original rank for ordering */
		if(group_id == 0){
			MPI_Comm_split(MPI_COMM_WORLD, group_id, myid_world, &MPI_COMM_GROUP);
		}else{
			MPI_Comm_split(MPI_COMM_WORLD, group_id, myid_world, &MPI_COMM_GROUP_OBSOLETE);
		}

		// barrier MPI_COMM_GROUP_OBSOLETE
		if(group_id != 0){

			MPI_Barrier(MPI_COMM_GROUP_OBSOLETE);

		}else{// MPI_COMM_GROUP

			/* get the size and rank in the new communicator */
			MPI_Comm_size(MPI_COMM_GROUP, &numprocs_group);
			MPI_Comm_rank(MPI_COMM_GROUP, &myid_group);

			// MPI_Cart
			// 2-d
			int ndims = 2;
			// each dimension size
			int dims[2];
			// each dimension is torus
			int periods[2] = {1,1};
			// do not reorder rank
			int reorder = 0;
			// coords of each element
			int coords[2];
			// cart size
			int cart_size = (int)sqrt( (double)numprocs_group );
			// neighbors
			int source_updown, destination_updown, source_leftright, destination_leftright;
			// each dimension size
			dims[0] = cart_size;
			dims[1] = numprocs_group/cart_size;

			/* print out original id, group id and dimensions in the group */
			printf("Numprocs_group: %i; Myid_world: %i; Myid_group: %i; Group: %i; dims[0]: %d; dims[1]: %d\n", numprocs_group, myid_world, myid_group, group_id, dims[0], dims[1]);
			printf("-------------------------------------\n");

			// cart communicator
			MPI_Comm MPI_COMM_CART;

			// create cart
			MPI_Cart_create(MPI_COMM_GROUP, ndims, dims, periods, reorder, &MPI_COMM_CART);

			//barrier useless processors out of the cart
			if( myid_group > (dims[0] * dims[1] - 1) ){

				printf("Useless Myid_group: %i\n", myid_group);
				printf("-------------------------------------\n");
				MPI_Barrier(MPI_COMM_GROUP);

			}else{

				MPI_Comm_size(MPI_COMM_CART, &numprocs);
				
				MPI_Comm_rank(MPI_COMM_CART, &myid);

				MPI_Cart_coords(MPI_COMM_CART, myid, ndims, coords);

				// Define the id of previous and next processor
				// Dimension 0, step is 1
				MPI_Cart_shift(MPI_COMM_CART, 0, 1, &source_updown, &destination_updown);
				// Dimension 1, step is 1
				MPI_Cart_shift(MPI_COMM_CART, 1, 1, &source_leftright, &destination_leftright);
				printf("Numprocs_cart: %i; Myid_cart: %i; source_leftright: %d; destination_leftright: %d\n", numprocs, myid, source_leftright, destination_leftright);
				printf("Numprocs_cart: %i; Myid_cart: %i; source_updown: %d; destination_updown: %d\n", numprocs, myid, source_updown, destination_updown);

				// Define edges and size of this processor updown and leftright
				int PROC_SIZE_UPDOWN = ( t/dims[0] - ( coords[0]<(t % dims[0])?0:1 ) + 1 ) * tile_size;
				int PROC_SIZE_LEFTRIGHT = ( t/dims[1] - ( coords[1]<(t % dims[1])?0:1 ) + 1 ) * tile_size;
				// print processor size for testing
				printf("Numprocs_cart: %i; Myid_cart: %i; PROC_SIZE_UPDOWN: %d; PROC_SIZE_LEFTRIGHT: %d\n", numprocs, myid, PROC_SIZE_UPDOWN, PROC_SIZE_LEFTRIGHT);
				printf("-------------------------------------\n");

				// Create grid
				// Each processor create its own blocks
				create_grid(&grid, PROC_SIZE_LEFTRIGHT, PROC_SIZE_UPDOWN);

				// Initialize grid
				board_init(&grid, PROC_SIZE_LEFTRIGHT, PROC_SIZE_UPDOWN);
				
				// Display grid
				// display_grid(grid, PROC_SIZE_LEFTRIGHT, PROC_SIZE_UPDOWN);

				// Four buffers for each processor
				int *bufferUp = malloc( sizeof(int) * PROC_SIZE_LEFTRIGHT );
				int *bufferDown = malloc( sizeof(int) * PROC_SIZE_LEFTRIGHT );
				int *bufferLeft = malloc( sizeof(int) * PROC_SIZE_UPDOWN );
				int *bufferRight = malloc( sizeof(int) * PROC_SIZE_UPDOWN );

				// array for storing left and right column for sending
				int *arrayLeft = malloc( sizeof(int) * PROC_SIZE_UPDOWN );
				int *arrayRight = malloc( sizeof(int) * PROC_SIZE_UPDOWN );

				// Computation starts
				while ( !finished && n_itrs < MAX_ITRS ){
					n_itrs++;
					
					int t, u, v;

					// transfer left and right boundary to array for sending
					for (t = 0; t < PROC_SIZE_UPDOWN; t++){
						arrayRight[t] = grid[t][PROC_SIZE_LEFTRIGHT-1];
						arrayLeft[t] = grid[t][0];
					}

					// 1. send maximum edge to bufferLeft of next processor, and get maximum edge of previous processor into its own bufferLeft
					MPI_Sendrecv(arrayRight, PROC_SIZE_UPDOWN, MPI_INT, destination_leftright, 1, 
						bufferLeft, PROC_SIZE_UPDOWN, MPI_INT, source_leftright, 1, MPI_COMM_CART, &status);

					// 2. send minimum edge to bufferRight of previous processor, and get minimum edge of next proccessor into its own bufferRight
					MPI_Sendrecv(arrayLeft, PROC_SIZE_UPDOWN, MPI_INT, source_leftright, 1, 
						bufferRight, PROC_SIZE_UPDOWN, MPI_INT, destination_leftright, 1, MPI_COMM_CART, &status);

					// 3. compare the minimum edge with bufferLeft and move red squares
					for (u =  0; u < PROC_SIZE_UPDOWN; u++){
						if (bufferLeft[u] == 1 && grid[u][0] == 0){
							grid[u][0] = 3;
							bufferLeft[u] = 4;
						}
					}

					// 4. move red squares in each processor
					/* red color movement */
					for (i =  0; i < PROC_SIZE_UPDOWN; i++){
						if (grid[i][0] == 1 && grid[i][1] == 0){
							grid[i][0] = 4;
							grid[i][1] = 3;
						}
						for (j = 1; j < PROC_SIZE_LEFTRIGHT; j++){
							if (grid[i][j] == 1 && (grid[i][(j+1)%PROC_SIZE_LEFTRIGHT] == 0)){
								grid[i][j] = 0;
								grid[i][(j+1)%PROC_SIZE_LEFTRIGHT] = 3;
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

					// 5. compare the maximum edge with bufferRight and move red squares
					for (v =  0; v < PROC_SIZE_UPDOWN; v++){
						if (bufferRight[v] == 0 && grid[v][PROC_SIZE_LEFTRIGHT-1] == 1){
							grid[v][PROC_SIZE_LEFTRIGHT-1] = 4;
							bufferLeft[v] = 3;
						}
					}

					// 6. send maximum edge to bufferUp of next processor, and get maximum edge of previous processor into its own bufferUp
					MPI_Sendrecv(grid[PROC_SIZE_UPDOWN-1], PROC_SIZE_LEFTRIGHT, MPI_INT, destination_updown, 1, 
						bufferUp, PROC_SIZE_LEFTRIGHT, MPI_INT, source_updown, 1, MPI_COMM_CART, &status);

					// 7. send minimum edge to bufferDown of previous processor, and get minimum edge of next proccessor into its own bufferDown
					MPI_Sendrecv(grid[0], PROC_SIZE_LEFTRIGHT, MPI_INT, source_updown, 1, 
						bufferDown, PROC_SIZE_LEFTRIGHT, MPI_INT, destination_updown, 1, MPI_COMM_CART, &status);

					// 8. compare the minimum edge with bufferUp and move blue squares
					for (k = 0; k < PROC_SIZE_LEFTRIGHT; k++)
					{
						if (bufferUp[k] == 2 && grid[0][k]==0){
							bufferUp[k] = 4;
							grid[0][k] = 3;
						}
					}

					// 9. move the blue squares of processor
					/* blue color movement */
					for (j = 0; j < PROC_SIZE_LEFTRIGHT; j++){
						if (grid[0][j] == 2 && grid[1][j] == 0){
							grid[0][j] = 4;
							grid[1][j] = 3;
						}
						// The maximum edge cannot move to its next rows otherwise array will get out of range. So the maximum cycle is 'PROC_SIZE - 1'
						for (i = 0 ; i < PROC_SIZE_UPDOWN - 1; i++){
							if (grid[i][j] == 2 && grid[(i+1)%PROC_SIZE_UPDOWN][j]==0){
								grid[i][j] = 0;
								// 这行的(i+1)%n应该为(i+1)%PROC_SIZE
								grid[(i+1)%PROC_SIZE_UPDOWN][j] = 3;
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

					// 10. compare the maximum edge with bufferDown and move blue squares
					for (l = 0; l < PROC_SIZE_LEFTRIGHT; l++)
					{
						if (bufferDown[l] == 0 && grid[PROC_SIZE_UPDOWN-1][l]==2){
							bufferDown[l] = 3;
							grid[PROC_SIZE_UPDOWN-1][l] = 4;
						}
					}

					// 11. check colored squares in each tile of processor
					int p, q, r, s;
					// loop in tiles in each row and colomn
					for (p = 0; p < ( PROC_SIZE_UPDOWN/tile_size ); p++){ // each column has p blocks
						for (q = 0; q < ( PROC_SIZE_LEFTRIGHT/tile_size ); q++){ // each row has q blocks

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
					
					// 12. Reduce all the 'finished' flag and do the logical OR operation.
					// Thus, if any processor has set 'finished' flag to TRUE, the result of logical OR will be TRUE.
					// Then, broadcast this 'finished' flag to each processors to terminate programs.
					MPI_Allreduce(&finished, &finished, 1, MPI_INT, MPI_LOR, MPI_COMM_CART);
					
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
				free(bufferLeft);
				free(bufferRight);
				free(arrayLeft);
				free(arrayRight);

				printf("Processor %d Program finished\n", myid);


				MPI_Barrier(MPI_COMM_GROUP);

			}// end else MPI_Barrier(MPI_COMM_GROUP);
			// printf("Myid_world: %d end else MPI_COMM_GROUP\n", myid_world);

		}// end else MPI_Barrier(MPI_COMM_GROUP_OBSOLETE);
		// printf("Myid_world: %d end else MPI_COMM_GROUP_OBSOLETE\n", myid_world);

	} // end parallel program
	// printf("Myid_world: %d end else parallel\n", myid_world);

	// MPI finalize
	MPI_Finalize();

	return 0;
}