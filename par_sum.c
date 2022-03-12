/*
 * par_sum.c
 *
 * CS 470 Project 1 (Pthreads)
 * Serial version
 *
 * @authors Teddy Pugh, Jacob Zornaik
 * Compile with --std=c99
 */

#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>

// aggregate variables
long sum = 0;
long odd = 0;
long min = INT_MAX;
long max = INT_MIN;
bool done = false;

pthread_mutex_t global_agg_m;
pthread_mutex_t task_queue_m;
pthread_cond_t cv = PTHREAD_COND_INITIALIZER;

bool done;

//function prototypes
void *worker(void *arg);
void update(long num);
void add_task(long num);
long remove_task();

// Representation of a task for the list.
typedef struct task {
    long val;
    struct task * next;
} node_t;

typedef struct list {
    node_t *head;
} queue_list;

queue_list task_queue;

void *worker(void *arg)
{
    //All threads skip this while condition
    while(task_queue.head != NULL || !done)
    {
        if(task_queue.head != NULL)
        {
            pthread_mutex_lock(&task_queue_m);
            long num = remove_task();
            printf("%ld\n", num);
            pthread_mutex_unlock(&task_queue_m);
            update(num);
        }
        else
        {
            pthread_cond_wait(&cv, &task_queue_m);
        }
    }
    return NULL;
}
/*
 * update global aggregate variables given a number
 */
void update(long num)
{
    // simulate computation
    sleep(num);

    // update aggregate variables
    pthread_mutex_lock(&global_agg_m);
    sum += num;
    //printf("%ld\n" , sum);
    pthread_mutex_unlock(&global_agg_m);
    if (num % 2 == 1) {
        pthread_mutex_lock(&global_agg_m);
        odd++;
        //printf("%ld\n" , odd);
        pthread_mutex_unlock(&global_agg_m);
    }
    if (num < min) {
        pthread_mutex_lock(&global_agg_m);
        min = num;
        pthread_mutex_unlock(&global_agg_m);
    }
    if (num > max) {
        pthread_mutex_lock(&global_agg_m);
        max = num;
        pthread_mutex_unlock(&global_agg_m);
    }
}

int main(int argc, char* argv[])
{
    // check and parse command line options
    if (argc != 3) {
        printf("Usage: sum <infile>\n");
        exit(EXIT_FAILURE);
    }
    char *fn = argv[1];
    int thread_count = atoi(argv[2]);
    
    // open input file
    FILE* fin = fopen(fn, "r");
    if (!fin) {
        printf("ERROR: could not open %s\n", fn);
        exit(EXIT_FAILURE);
    }
    //initialize work queue and sync variables
    pthread_t workers[thread_count];
    //spawn worker threads based on command line args, for now goal
    // is to get it working with just 1 worker thread, then add functionality
    // for multiple threads (change pthread_t workers to an array with malloc).
    for(int i = 0; i < thread_count; i++)
    {
        pthread_create(&workers[i], NULL, worker, NULL);
    }
    //for each (action,num) pair in input:
    //    //if action == p
    //        add num to work queue, wake an idle worker thread
    //    else if action == w
    //         wait num seconds
    //
    //done = true
    //exhaust work queue and wait for workers to finish
    //

    
    // load numbers and add them to the queue
    char action;
    long num;
    while (fscanf(fin, "%c %ld\n", &action, &num) == 2) {
        if (action == 'p') {            // process
            add_task(num);
            pthread_cond_signal(&cv);
        } else if (action == 'w') {     // wait
            sleep(num);
        } else {
            printf("ERROR: Unrecognized action: '%c'\n", action);
            exit(EXIT_FAILURE);
        }
    }

    done = true;
    fclose(fin);
    for(int i = 0; i < thread_count; i++)
    {
        pthread_join(workers[i], NULL);
    }
    // print results
    printf("%ld %ld %ld %ld\n", sum, odd, min, max);

    // clean up and return
    pthread_cond_destroy(&cv);
    pthread_mutex_destroy(&global_agg_m);
    pthread_mutex_destroy(&task_queue_m);
    return (EXIT_SUCCESS);
}

void add_task(long num)
{
    node_t new_task;
    new_task.val = num;
    // if there are no elements in the list, add the task
    if(task_queue.head == NULL)
    {
        task_queue.head = &new_task;
    }
    // otherwise, add it to the list of tasks
    else
    {
        node_t* curr = task_queue.head;
        while(curr->next != NULL)
        {
            curr = curr->next;
        }
        curr->next = &new_task;
    }
}

long remove_task()
{
    long value;
    if(task_queue.head == NULL)
    {
        printf("ERROR, not supposed to be here.\n");
        exit(EXIT_FAILURE);
    }
    else if(task_queue.head->next == NULL)
    {
        value = task_queue.head->val;
        task_queue.head = NULL;
    }
    else
    {
        value = task_queue.head->val;
        task_queue.head = task_queue.head->next;
    }

    return value;
}

