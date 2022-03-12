/*
 * mergesort.c
 *
 * CS 470 Project 2 (MPI)
 * Original serial version.
 *
 * Name(s): Teddy Pugh and Jacob Zorniak
 */

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <mpi.h>

// histogram bins
#define BINS 10

// maximum random number
#define RMAX 100

// enable debug output
//#define DEBUG

// timing macros (must first declare "struct timeval tv")
#define START_TIMER(start) start = MPI_Wtime();
#define STOP_TIMER(end) end = MPI_Wtime();
#define GET_TIMER(start,end) (end-start)

// "count_t" used for number counts that could become quite high
typedef unsigned long count_t;

int *nums;            // random numbers
int *local_nums;      // local copy of nums for each process
count_t  global_n;    // global "nums" count
count_t  local_n;     // local "local_nums" count
count_t  shift_n;     // global left shift offset
count_t *hist;        // histogram (counts of "nums" in bins)

// varialbes for determining rank or total number of processes
int mpi_rank;
int mpi_size;

/*
 *  * Helper function for providing a log base 2 operation. 
 *   * Given value must be a power of 2.
 *    */  
int lg2(int a)
{   
    int result = 0;
    while(a > 1)
    {
        a = a/2;
        result++;
    }
    return result;
}

/*
 *  * Helper function for providing a exponentiation operation.
 *   */   
int Pow(int a, int b)
{
    int result = 1;
    for(int i = 1; i <= b; i++)
    {
        result = result * a;
    }
    return result;
}
/*
 * Helper function for providing a modulus operation.
 */
int mod (int a, int b)
{
    int ret = a%b;
    if(ret < 0)
    {
        ret += b;
    }
    return ret;
}

/*
 * Parse and handle command-line parameters. Returns true if parameters were
 * valid; false if not.
 */
bool parse_command_line(int argc, char *argv[])
{
    // read command-line parameters
    if (argc != 3) {
        printf("Usage: %s <n> <shift>\n", argv[0]);
        return false;
    } else {
        global_n = strtol(argv[1], NULL, 10);
        shift_n  = strtol(argv[2], NULL, 10);
    }

    // check shift offset
    if (shift_n > global_n) {
        printf("ERROR: shift offset cannot be greater than N\n");
        return false;
    }

    return true;
}

/*
 * Allocate and initialize number array and histogram.
 */
void initialize_data_structures()
{
    // initialize local data structures
    nums = (int*)calloc(global_n, sizeof(int));
    if (nums == NULL) {
        fprintf(stderr, "Out of memory!\n");
        exit(EXIT_FAILURE);
    }
    hist = (count_t*)calloc(BINS, sizeof(count_t));
    if (hist == NULL) {
        fprintf(stderr, "Out of memory!\n");
        exit(EXIT_FAILURE);
    }
}

/*
 * Compares two ints. Suitable for calls to standard qsort routine.
 */
int cmp(const void* a, const void* b)
{
    return *(int*)a - *(int*)b;
}

/*
 * Print contents of an int list.
 */
void print_nums(int *a, count_t n)
{
    for (count_t i = 0; i < n; i++) {
        printf("%d ", a[i]);
    }
    printf("\n");
}

/*
 * Print contents of a count list (i.e., histogram).
 */
void print_counts(count_t *a, count_t n)
{
    for (count_t i = 0; i < n; i++) {
        printf("%lu ", a[i]);
    }
    printf("\n");
}

/*
 * Merge two sorted lists ("left" and "right) into "dest" using temp storage.
 */
void merge(int left[], count_t lsize, int right[], count_t rsize, int dest[])
{
    count_t dsize = lsize + rsize;
    int *tmp = (int*)malloc(sizeof(int) * dsize);
    if (tmp == NULL) {
        fprintf(stderr, "Out of memory!\n");
        exit(EXIT_FAILURE);
    }
    count_t l = 0, r = 0;
    for (count_t ti = 0; ti < dsize; ti++) {
        if (l < lsize && (left[l] <= right[r] || r >= rsize)) {
            tmp[ti] = left[l++];
        } else {
            tmp[ti] = right[r++];
        }
    }
    memcpy(dest, tmp, dsize*sizeof(int));
    free(tmp);
}

/*
 * Generate random integers for "nums".
 * 
 * Analysis:
 * In order to parallelize this function using MPI calls, we simply filled the num[] with random numbers up to 100. 
 * Then, we used an MPI_Scatter call to distribute the random numbers evenly to each process. This subroutine functions correctly
 * when run with 64 processes and 256 nums. However, it does not demostrate strong nor weak scaling. As the number of processes increases,
 * for both a static and dynamic number of random numbers, the run-time increases. The same trends can be observed when running the subroutine
 * with the same number of processe and only changing the amount of random numbers. 
 *
 * The reason for this type of performance lies in the fact that the subroutine is essentially running the program in serial and then
 * calling MPI_Scatter to distribute the random numbers across the processes. Thus, this subroutine is more communication-bound as the 
 * same computation is being done to calculate the same number of random values but since we increase the amount of processes, we
 * increase the amount of communcation that must be done by MPI_Scatter which may even have to be across nodes.
 *
 * Graph: https://docs.google.com/spreadsheets/d/18GJocgaWcoZSSXl5I5JyTfpNLQrsHDU4OrSGAUP5mck/edit?usp=sharing
 *
 * +---------------------+------------------+----------+
 * | Number of Processes | Number of values | Run-time |
 * +_____________________+__________________+__________+
 * |     serial          |    130000000     |  1.7348  |
 * +_____________________+__________________+__________+
 * |      4              |    130000000     |  2.9745  |
 * +_____________________+__________________+__________+
 * |      8              |    130000000     |  3.2569  |
 * +_____________________+__________________+__________+
 * |      16             |    130000000     |  5.4111  |
 * +_____________________+__________________+__________+
 * |      32             |    130000000     |  6.4867  |
 * +_____________________+__________________+__________+
 * |      64             |    130000000     |  7.0242  |
 * +_____________________+__________________+__________+
 * +_____________________+__________________+__________+
 * |      4              |    4000000       |  0.1124  |
 * +_____________________+__________________+__________+
 * |      8              |    8000000       |  0.2057  |
 * +_____________________+__________________+__________+
 * |      16             |    16000000      |  0.6722  |
 * +_____________________+__________________+__________+
 * |      32             |    32000000      |  1.5765  |
 * +_____________________+__________________+__________+
 * |      64             |    64000000      |  3.4621  |
 * +_____________________+__________________+__________+
 * +_____________________+__________________+__________+
 * |      4              |    4000000       |  0.1124  |
 * +_____________________+__________________+__________+
 * |      4              |    8000000       |  0.1778  |
 * +_____________________+__________________+__________+
 * |      4              |    16000000      |  0.3992  |
 * +_____________________+__________________+__________+
 * |      4              |    32000000      |  0.8536  |
 * +_____________________+__________________+__________+
 * |      4              |    64000000      |  1.4850  |
 * +_____________________+__________________+__________+
 *
 */
void randomize()
{
    // Generate all the random values on process 0
    if(mpi_rank == 0)
    {
        srand(42);
        for (count_t i = 0; i < global_n; i++) 
        {
            nums[i] = rand() % RMAX;
        }
    }
    // distribute the random values evenly to each process
    MPI_Scatter(nums,       local_n, MPI_INT, 
                local_nums, local_n, MPI_INT, 0, MPI_COMM_WORLD);
}

/*
 * Calculate histogram based on contents of "nums".
 *
 * Analysis: To parallelize this subroutine, we split up the work across processes to create their own histogram
 * and then we perform a reduction to aggregate the results to process 0. The subroutine functions correctly when run with 64 processes
 * and 256 random numbers. This subroutine demonstrates strong scaling up to 64 processes and 130 million random numbers and actually has 
 * a linear speedup when using 8 or more processes. As can be observed in the table below, when the number of processes doubles (after 8 processes),
 * the run-time is cut in half. This subroutine also demonstrates weak scaling as doubling the number of random numbers and the number of processes results 
 * in the run-time staying roughly the same. Also, doubling the problem size while not changing the number of processes results in double the run-time which
 * further proves linear speedup.
 *
 * The reason this subroutine manages to perform with linear speedup is due to the efficiency of parallelization in splitting up the problem
 * for each process and communicating using MPI_Reduce in an efficient manner. This subroutine is more computation bound as the communication
 * works very efficiently and is simply held back by the amount of time it takes for each process to complete its computation.
 *
 * Graph: https://docs.google.com/spreadsheets/d/1VGQmkZGhd7kgrFB-wbiw8yaiB2cykiJvsBEqeioFd_g/edit?usp=sharing 
 *
 * +---------------------+------------------+----------+
 * | Number of Processes | Number of values | Run-time |
 * +_____________________+__________________+__________+
 * |     serial          |    130000000     |  0.4205  |
 * +_____________________+__________________+__________+
 * |      4              |    130000000     |  0.0645  |
 * +_____________________+__________________+__________+
 * |      8              |    130000000     |  0.0715  |
 * +_____________________+__________________+__________+
 * |      16             |    130000000     |  0.0352  |
 * +_____________________+__________________+__________+
 * |      32             |    130000000     |  0.0170  |
 * +_____________________+__________________+__________+
 * |      64             |    130000000     |  0.0089  |
 * +_____________________+__________________+__________+
 * +_____________________+__________________+__________+
 * |      4              |    4000000       |  0.0037  |
 * +_____________________+__________________+__________+
 * |      8              |    8000000       |  0.0044  |
 * +_____________________+__________________+__________+
 * |      16             |    10000000      |  0.0043  |
 * +_____________________+__________________+__________+
 * |      32             |    30000000      |  0.0045  |
 * +_____________________+__________________+__________+
 * |      64             |    60000000      |  0.0045  |
 * +_____________________+__________________+__________+
 * +_____________________+__________________+__________+
 * |      4              |    4000000       |  0.0037  |
 * +_____________________+__________________+__________+
 * |      4              |    8000000       |  0.0074  |
 * +_____________________+__________________+__________+
 * |      4              |    16000000      |  0.0086  |
 * +_____________________+__________________+__________+
 * |      4              |    32000000      |  0.0159  |
 * +_____________________+__________________+__________+
 * |      4              |    64000000      |  0.0317  |
 * +_____________________+__________________+__________+
 *
 */
void histogram()
{
    // CALCULATE the local histogram
    for (count_t i = 0; i < local_n; i++) {
        hist[local_nums[i] % BINS]++;
    }
    
    // AGGREGATE results to process 0
    count_t* temp = (count_t*)malloc(BINS*sizeof(count_t)); // temporary variable for storing the reduce
    MPI_Reduce(hist, temp, BINS, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD); 
    if(mpi_rank == 0) // copy the results from temp back to hist in process 0
    {
        hist = temp;
    }
}

/*
 * Shift "nums" left by the given number of slots. Anything shifted off the left
 * side should rotate around to the end, so no numbers are lost.
 *
 * Analysis: This subroutine functions correctly when run with 64 processes and 256 random numbers. It does not demonstrate strong nor weak scaling.
 * As the number of processes increases, the run-time increases. Also, as the number of processes increases and the problem size increases, the run-time increases. 
 * There are a few jumps in run-time that can be observed. For example, between 8 to 16 processes, 16 to 32 processes, and 32 to 64 processes. 
 * This may be caused due to an increased number of nodes being used for the processes and the increase in communication time between
 * these nodes. This would imply that the subroutine is more communication-bound.
 *
 * Changing the shift amount while maintaining the same number of random numbers and processes resulted in roughly the same run-time for every shift amount.
 * This subroutine functions in this way because it is slowed down by the send and receive funtions which may happen across nodes in some cases.
 *
 * Graph: https://docs.google.com/spreadsheets/d/1TM0nuFRrZp6G1mYZH1d4WBm_slBsK8wE25LdGvJaH0Y/edit?usp=sharing
 *
 * +---------------------+-------------------------+----------+
 * | Number of Processes | Number of values, shift | Run-time |
 * +_____________________+_________________________+__________+
 * |    serial           |    130000000, 2         |  1.3628  |
 * +_____________________+_________________________+__________+
 * |      4              |    130000000, 2         |  0.1107  |
 * +_____________________+_________________________+__________+
 * |      8              |    130000000, 2         |  0.1187  |
 * +_____________________+_________________________+__________+
 * |      16             |    130000000, 2         |  2.2346  |
 * +_____________________+_________________________+__________+
 * |      32             |    130000000, 2         |  3.3278  |
 * +_____________________+_________________________+__________+
 * |      64             |    130000000, 2         |  3.8757  |
 * +_____________________+_________________________+__________+
 * +_____________________+_________________________+__________+
 * |      4              |     4000000, 2          |  0.0061  |
 * +_____________________+_________________________+__________+
 * |      8              |     8000000, 2          |  0.0094  |
 * +_____________________+_________________________+__________+
 * |      16             |     16000000, 2         |  0.2790  |
 * +_____________________+_________________________+__________+
 * |      32             |     32000000, 2         |  0.8219  |
 * +_____________________+_________________________+__________+
 * |      64             |     64000000, 2         |  1.9100  |
 * +_____________________+_________________________+__________+
 * +_____________________+_________________________+__________+
 * |      64             |    1000000, 4           |  0.0317  |
 * +_____________________+_________________________+__________+
 * |      64             |    1000000, 8           |  0.0316  |
 * +_____________________+_________________________+__________+
 * |      64             |    1000000, 16          |  0.0317  |
 * +_____________________+_________________________+__________+
 * |      64             |    1000000, 32          |  0.0318  |
 * +_____________________+_________________________+__________+
 * |      64             |    1000000, 64          |  0.0308  |
 * +_____________________+_________________________+__________+
 *
 */
void shift_left()
{
    // preserve first shift_n values
    int *tmp = (int*)malloc(sizeof(int) * local_n);
    if (tmp == NULL) {
        fprintf(stderr, "Out of memory!\n");
        exit(EXIT_FAILURE);
    }
    // copy the local_nums to the temporary array.
    memcpy(tmp, local_nums, sizeof(int)*local_n);
    // shift the tmp array by the shift amount
    MPI_Sendrecv(tmp                          , shift_n, MPI_INT, mod((mpi_rank-1), mpi_size), 0,
                 &local_nums[local_n-shift_n], shift_n, MPI_INT, MPI_ANY_SOURCE,              0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // copy the non shifted values from the tmp array back to local
    MPI_Sendrecv(&tmp[shift_n],  mpi_size-shift_n, MPI_INT, mpi_rank,       0,
                 &local_nums[0], mpi_size-shift_n, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // gather the values back to the root process's nums
    MPI_Gather(local_nums, local_n, MPI_INT,
               nums,       local_n, MPI_INT, 0, MPI_COMM_WORLD);
    free(tmp);
}

/*
 * Merge sort helper (shouldn't be necessary in parallel version).
 */
void merge_sort_helper(int *start, count_t len)
{
    if (len < 100) {
        qsort(start, len, sizeof(int), cmp);
    } else {
        count_t mid = len/2;
        merge_sort_helper(start, mid);
        merge_sort_helper(start+mid, len-mid);
        merge(start, mid, start+mid, len-mid, start);
    }
}

/*
 * Sort "nums" using the mergesort algorithm.
 *
 * Analysis: To parallelize this subroutine we used a MPI_Scatter call to distribute the nums to each process. Then each process performs a quick sort on their local_nums.
 * After that, we calculate the height of the tree used for communication (sending and receiving). We use a for loop from i = 0 to i < height. The nodes that are divisible
 * by 2^(i + 1) will be receiving local_nums from a sending node and will merge the two arrays into one. The nodes that are divisible by 2^i will be sending their local_nums array.
 * In the end we are left with the merge sorted array that we copy from local_nums to nums using a temp pointer so we can free nums after the copy is complete.
 *
 * This subroutine functions correctly when run with 64 processes and 256 random numbers. There is clearly a speedup in run-time from the serial version of the subroutine.
 * We can also observe strong scaling (linear speedup) and weak scaling when the subroutine is run from 4 to 8 processes. However, after 8 processes the run-time increases as the number
 * of processes doubles. This could be caused by the same problem we observed in the shift_left() subroutine where once the number of processes passes 8, more nodes must be used
 * which slows down communication between the send and receive functions. Thus, this subroutine is mainly communication-bound.
 *
 * Graph: https://docs.google.com/spreadsheets/d/1CzwQCL7i3n-7gvm1uyZsk8QiCOZtBit039oCcBRP9iE/edit?usp=sharing 
 *
 * +---------------------+------------------+----------+
 * | Number of Processes | Number of values | Run-time |
 * +_____________________+__________________+__________+
 * |     serial          |     130000000    |  12.2295 |
 * +_____________________+__________________+__________+
 * |      4              |     130000000    |  4.1241  |
 * +_____________________+__________________+__________+
 * |      8              |     130000000    |  2.6252  |
 * +_____________________+__________________+__________+
 * |      16             |     130000000    |  6.0022  |
 * +_____________________+__________________+__________+
 * |      32             |     130000000    |  7.6749  |
 * +_____________________+__________________+__________+
 * |      64             |     130000000    |  8.5490  |
 * +_____________________+__________________+__________+
 * +_____________________+__________________+__________+
 * |      4              |     4000000      |  0.1157  |
 * +_____________________+__________________+__________+
 * |      8              |     8000000      |  0.1510  |
 * +_____________________+__________________+__________+
 * |      16             |     16000000     |  0.7315  |
 * +_____________________+__________________+__________+
 * |      32             |     32000000     |  1.8969  |
 * +_____________________+__________________+__________+
 * |      64             |     64000000     |  4.2092  |
 * +_____________________+__________________+__________+
 */
void merge_sort()
{
    //distribute the nums to each process
    MPI_Scatter(nums,       local_n, MPI_INT,
        local_nums, local_n, MPI_INT, 0, MPI_COMM_WORLD);
    // sort each procceses local_nums using qsort
    qsort(local_nums, local_n, sizeof(int), cmp);
    int height = lg2(mpi_size);
    for(int i = 0; i < height; i++)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        if( (mpi_rank % (Pow(2,(i+1)))) == 0)
        {
            int *tmp = (int*)malloc(sizeof(int) * global_n);
            MPI_Recv(tmp, global_n, MPI_INT, mpi_rank + Pow(2,i), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            merge(local_nums, local_n * Pow(2,i), tmp, local_n * Pow(2,i), local_nums); 
        }
        else if(mpi_rank % ((Pow(2,i))) == 0)
        {
            MPI_Send(local_nums, local_n * Pow(2, i), MPI_INT, mpi_rank - Pow(2,i), 0, MPI_COMM_WORLD); 
        }
        // if it does not enter, then the value has already been merged.      
        MPI_Barrier(MPI_COMM_WORLD);
    }
    if(mpi_rank == 0)
    {
        //create temporary pointer for freeing old nums
        int* tmp;
        tmp = nums;
        nums = local_nums;
        free(tmp);
    }
}

int main(int argc, char *argv[])
{

    if (!parse_command_line(argc, argv)) {
        exit(EXIT_FAILURE);
    }
    
    //prepare MPI related variables & functions
    MPI_Init(&argc, &argv);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    initialize_data_structures();
    
    // instantiate local nums to be used by the non root processes
    local_n = global_n/mpi_size;
    local_nums = (int*)(malloc(sizeof(int)*global_n));

    double start_time;
    double end_time;
    double rand_time;
    double hist_time;
    double shft_time;
    double sort_time;

    START_TIMER(start_time);
    randomize();
    MPI_Barrier(MPI_COMM_WORLD);
    STOP_TIMER(end_time);
    rand_time = GET_TIMER(start_time, end_time);

#   ifdef DEBUG
    if(mpi_rank == 0) { printf("global orig list: "); print_nums(nums, global_n); }
#   endif

    // compute histogram
    START_TIMER(start_time);
    histogram();
    MPI_Barrier(MPI_COMM_WORLD);
    STOP_TIMER(end_time);
    hist_time = GET_TIMER(start_time, end_time);

    // print histogram
    if(mpi_rank == 0) 
    {
        printf("GLOBAL hist: "); print_counts(hist, BINS);
    }

    // perform left shift
    START_TIMER(start_time);
    shift_left();
    MPI_Barrier(MPI_COMM_WORLD);
    STOP_TIMER(end_time);
    shft_time = GET_TIMER(start_time, end_time);

#   ifdef DEBUG
    if(mpi_rank == 0) { printf("global shft list: "); print_nums(nums, global_n); }
#   endif

    // perform merge sort
    START_TIMER(start_time);
    merge_sort();
    MPI_Barrier(MPI_COMM_WORLD);
    STOP_TIMER(end_time);
    sort_time = GET_TIMER(start_time, end_time);

    // print global results
#   ifdef DEBUG
    if(mpi_rank == 0) { printf("GLOBAL list: "); print_nums(nums, global_n); }
#   endif
    
    if(mpi_rank == 0)
    {
        printf("RAND: %.4f  HIST: %.4f  SHFT: %.4f  SORT: %.4f\n",
                rand_time, hist_time, shft_time, sort_time);
    }
    
    // clean up and exit
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    free(local_nums); //nums was freed in merge
    free(hist);
    return EXIT_SUCCESS;
}
