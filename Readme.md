#COMP5426 PARALLEL AND DISTRIBUTED COMPUTING

##Assignment
---

###Assignment1 The Red/Blue computation (MPI)

The Red/Blue computation simulates two interactive flows: an n by n board is initialized so cells have one of three colors: red, white, and blue, where white is empty, red moves right, and blue moves down. (The board may be initialized with 1/3 cells in read, 1/3 in white and 1/3 in blue and colors should be interleaved and spread across the board. You need to write a separate function board_init to initialize the board.) Colors wraparound to the opposite side when reaching the edge. In the first half step of an iteration, any red color can move right one cell if the cell to the right is unoccupied (white); on the second half step, any blue color can move down one cell if the cell below it is unoccupied (white); the case where red vacates a cell (first half) and blue moves into it (second half) is okay. Viewing the board as overlaid with t by t tiles (t divides n), the computation terminates if any tile’s colored squares are more than c% one color (blue or red).

---

###Assignment2 Sieve of Eratosthenes (Pthreads)

Using Pthreads to write a program to find all prime numbers within a range from 2 to n using the Sieve of Eratosthenes.

Considering the load balance of threads. Each thread has access to part of the global array and will do sieve algorithm from 2 to sqrt(thread_end). Do not use the partition based on checking base number since it is not load balancing.