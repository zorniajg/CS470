/**
 * CURRENTLY DOES NOT RUN HOW IT IS SUPPOSED TO (HAS BUGS AND CODE WILL BE CHANGED SIGNIFICANTLY FOR FINAL DELIVERABLE).
 *
 * Parallel implementation of Dijkstra's shortest path algorithm with timing mechanisms.
 * This is not a new way to implement the alogrithm in a parallel manner, 
 * so the psuedocode and much of the code was sourced from these two sources
 * ( https://cse.buffalo.edu/faculty/miller/Courses/CSE633/Ye-Fall-2012-CSE633.pdf ,
 *   https://github.com/Lehmannhen/MPI-Dijkstra/commit/5991fb8c938c44e5b53da7210a49b28633d11989#diff-1e47380c4e074b2ea5da122e1d512c8dR429)
 *
 **/

#include<stdio.h>
#include<stdlib.h>
#include<sys/time.h>
#include<mpi.h>
#define INFINITY 1000000
#define MAX 10

#define START_TIMER(start) start = MPI_Wtime();
#define STOP_TIMER(end) end = MPI_Wtime();
#define GET_TIMER(start,end) (end-start)

void dijkstra(int *local_matrix, int *local_distance, int *local_pred, int local_n, int n);
MPI_Datatype build_blk_col_type(int n, int loc_n);
int find_min_dist(int loc_dist[], int loc_known[], int loc_n);
void print_dists(int global_dist[], int n);
void print_paths(int global_pred[], int n);

int *local_matrix, *local_distance, *local_pred;
int *global_distance, *global_pred;
int global_n, local_n;

int mpi_rank;
int mpi_size;

int main(int argc, char *argv[])
{
        int G[MAX][MAX],i,j,n;
        MPI_Datatype blk_col_mpi_t;      

        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
        
        if(mpi_rank == 0)
        {   
            printf("Enter no. of vertices:");
            scanf("%d",&n);
        }
        //allocate and intialize necessary structures and variables
       local_n = global_n/mpi_size;
       local_matrix = (int*)(malloc(sizeof(int) * local_n * n));
       local_distance = (int*)(malloc(sizeof(int) * local_n));
       local_pred = (int*)(malloc(sizeof(int) * local_n));
       blk_col_mpi_t = build_blk_col_type(n, local_n);

       if(mpi_rank == 0)
       {
           global_distance = (int*)(malloc(sizeof(int) * n));
           global_pred = (int*)(malloc(sizeof(int) * n));
       }
       
       //Get the matrix and scatter block columns to each process
       printf("\nEnter the adjacency matrix:\n");

            for(i=0;i<n;i++)
                for(j=0;j<n;j++)
                    scanf("%d",&G[i][j]);

       MPI_Scatter(G, 1, blk_col_mpi_t, local_matrix, n * local_n, MPI_INT, 0, MPI_COMM_WORLD);
       
        //Initialize timing variables
        double start_time;
        double end_time;
        double dijkstra_time;
        
       //Get the run-time of the parallel Dijkstra's algorithm 
        START_TIMER(start_time);
        dijkstra(local_matrix, local_distance, local_pred, local_n, n);

        MPI_Barrier(MPI_COMM_WORLD);
        STOP_TIMER(end_time);
        dijkstra_time = GET_TIMER(start_time, end_time);
        printf("s: %.4f e: %.4f d: %.4f \n", start_time, end_time, dijkstra_time);

        //print the path and distance of each node
        if(mpi_rank == 0)
        {
            //print_dists(global_distance, n);
            //print_paths(global_pred, n);
            printf("Run-time: %.4f\n", dijkstra_time);
            free(global_distance);
            free(global_pred);
        }

        free(local_matrix);
        free(local_pred);
        free(local_distance);
        MPI_Type_free(&blk_col_mpi_t);
        MPI_Finalize();
        return EXIT_SUCCESS;
}

void dijkstra(int *local_matrix, int *local_distance, int *local_pred, int local_n, int n)
{
        
        int *local_visited;
        int i, local_v, local_u, global_u, new_distance, distance_global_u;
        int my_min[2];
        int global_min[2];
       
        local_visited = (int*)(malloc(sizeof(int) * local_n)); 
        //initialize local data structures
        if(mpi_rank == 0)
            local_visited[0] = 1;
        else
            local_visited[0] = 0;

        for(local_v = 1; local_v < local_n; local_v++)
            local_visited[local_v] = 0;

        for(local_v = 0; local_v < local_n; local_v++)
        {
            local_distance[local_v] = local_matrix[0 * local_n + local_v];
            local_pred[local_v] = 0;
        }

        //Find min distance for local vertices
        for(i = 0; i < n - 1; i++)
        {
            local_u = find_min_dist(local_distance, local_visited, local_n);
            
            if(local_u != -1)
            {
                my_min[0] = local_distance[local_u];
                my_min[1] = local_u + mpi_rank * local_n;
            }
            else
            {
                my_min[0] = INFINITY;
                my_min[1] = -1;
            }
            //Get the min distance found by the processes and store that distance and the global vertex in the global_min
            MPI_Allreduce(my_min, global_min, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD);

            distance_global_u = global_min[0];
            global_u = global_min[1];

            if(global_u == -1)
                break;
            //Check if global u belongs to current process, and if so update local_visited
            if((global_u / local_n) == mpi_rank)
            {
                local_u = global_u % local_n;
                local_visited[local_u] = 1;
            }
            //Update the distances from the source vertex to local_v for each local vertex and check if an unvisited vertex's
            //distance from the source to global_u + distance from global_u to local V is smaller than the distance from the source to local v.
            for(local_v = 0; local_v < local_n; local_v++)
            {
                if(!local_visited[local_v])
                {
                    new_distance = distance_global_u + local_matrix[global_u * local_n + local_v];
                    if(new_distance < local_distance[local_v])
                    {
                        local_distance[local_v] = new_distance;
                        local_pred[local_v] = global_u;
                    }
                }
            }
       }
            free(local_visited);
}

/*---------------------------------------------------------------------
 * Function:  Build_blk_col_type
 * Purpose:   Build an MPI_Datatype that represents a block column of
 *            a matrix
 * In args:   n:  number of rows in the matrix and the block column
 *            loc_n = n/p:  number cols in the block column
 * Ret val:   blk_col_mpi_t:  MPI_Datatype that represents a block
 *            column
 */
MPI_Datatype build_blk_col_type(int n, int loc_n) {
    MPI_Aint lb, extent;
    MPI_Datatype block_mpi_t;
    MPI_Datatype first_bc_mpi_t;
    MPI_Datatype blk_col_mpi_t;

    MPI_Type_contiguous(loc_n, MPI_INT, &block_mpi_t);
    MPI_Type_get_extent(block_mpi_t, &lb, &extent);

    /* MPI_Type_vector(numblocks, elts_per_block, stride, oldtype, *newtype) */
    MPI_Type_vector(n, loc_n, n, MPI_INT, &first_bc_mpi_t);

    /* This call is needed to get the right extent of the new datatype */
    MPI_Type_create_resized(first_bc_mpi_t, lb, extent, &blk_col_mpi_t);
    MPI_Type_commit(&blk_col_mpi_t);
    MPI_Type_free(&block_mpi_t);
    MPI_Type_free(&first_bc_mpi_t);

    return blk_col_mpi_t;
}

/*-------------------------------------------------------------------
 * Function:   Find_min_dist
 * Purpose:    find the minimum local distance from the source to the
 *             assigned vertices of the process that calls the method
 *
 *
 * In args:    loc_dist:  array with distances from source 0
 *             loc_known: array with values 1 if the vertex has been visited
 *                        0 if not
 *             loc_n:     local number of vertices
 *
 * Return val: loc_u: the vertex with the smallest value in loc_dist,
 *                    -1 if all vertices are already known
 *
 * Note:       loc_u = -1 is not supposed to be used when this function returns
 *
 */
int find_min_dist(int loc_dist[], int loc_known[], int loc_n) {
    int loc_u = -1, loc_v;
    int shortest_dist = INFINITY;
    for (loc_v = 0; loc_v < loc_n; loc_v++) {
        if (!loc_known[loc_v]) {
            if (loc_dist[loc_v] < shortest_dist) {
                shortest_dist = loc_dist[loc_v];
                loc_u = loc_v;
            }
        }
    }
    return loc_u;
}

/*-------------------------------------------------------------------
 *  * Function:    Print_dists
 *   * Purpose:     Print the length of the shortest path from 0 to each
 *    *              vertex
 *     * In args:     n:  the number of vertices
 *      *              dist:  distances from 0 to each vertex v:  dist[v]
 *       *                 is the length of the shortest path 0->v
 *        */
void print_dists(int global_dist[], int n) {
    int v;
    printf("  v    dist 0->v\n");
    printf("----   ---------\n");
    for (v = 1; v < n; v++) {
        if (global_dist[v] == INFINITY) {
            printf("%3d       %5s\n", v, "inf");
        }
        else
            printf("%3d       %4d\n", v, global_dist[v]);
        }
    printf("\n");
}
/*-------------------------------------------------------------------
 *  * Function:    Print_paths
 *   * Purpose:     Print the shortest path from 0 to each vertex
 *    * In args:     n:  the number of vertices
 *     *              pred:  list of predecessors:  pred[v] = u if
 *      *                 u precedes v on the shortest path 0->v
 *       */
void print_paths(int global_pred[], int n) {
    int v, w, *path, count, i;
    path =  malloc(n * sizeof(int));
    printf("  v     Path 0->v\n");
    printf("----    ---------\n");
    for (v = 1; v < n; v++) {
        printf("%3d:    ", v);
        count = 0;
        w = v;
        while (w != 0) {
            path[count] = w;
            count++;
            w = global_pred[w];
        }
        printf("0 ");
        for (i = count-1; i >= 0; i--)
            printf("%d ", path[i]);
        printf("\n");
    }
    free(path);
}
