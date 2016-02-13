/*
  Original Author: Joshua Stough, Washington and Lee University This
  code was initially obtained from
  http://sc12.supercomputing.org/hpceducator/PythonForParallelism/codes/parallelQuicksort.c.
  Later, Arnaud Legrand mainly made a few cosmetic changes so that it
  is easier to use for performance evaluation purposes.

  This code quicksorts a random list of size given by the argument
  (default 1M) and times both sequential quicksort and parallel (using
  Pthreads).
*/


#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>

#include <algorithm>    // std::shuffle


#define DNUM 1000000
#define THREAD_LEVEL 3

#define REP 3


//for sequential and parallel implementation
void swap(T lyst[], int i, int j);
int partition(T lyst[], int lo, int hi);
void quicksortHelper(T lyst[], int lo, int hi);
void quicksort(T lyst[], int size);
int isSorted(T lyst[], int size);

//for parallel implementation
void parallelQuicksort(T lyst[], int size, int tlevel);
void *parallelQuicksortHelper(void *threadarg) __attribute__ ((noreturn));

struct thread_data {
  T *lyst;
  int low;
  int high;
  int level;
};
//thread_data should be thread-safe, since while lyst is
//shared, [low, high] will not overlap among threads.

//for the builtin qsort, for fun:
int compare_doubles(const void *a, const void *b);


struct timer1 {

  struct timeval _start, _end;
  void start() {
    gettimeofday(&_start, NULL);
  }
  void end() {
    gettimeofday(&_end, NULL);
  }
  double get_elapsed_time() {
    return ((_end.tv_sec * 1000000 + _end.tv_usec)
            - (_start.tv_sec * 1000000 + _start.tv_usec)) / 1000000.0;
  }
};

struct timer2 {

  #define BILLION 1E9
  // CLOCK_PROCESS_CPUTIME_ID, CLOCK_REALTIME

  struct timespec _start, _end;
  void start() {
    clock_gettime(CLOCK_REALTIME, &_start);
  }
  void end() {
    clock_gettime(CLOCK_REALTIME, &_end);
  }
  double get_elapsed_time() {
    return (_end.tv_sec - _start.tv_sec)
    + (_end.tv_nsec - _start.tv_nsec) / BILLION;
  }
};


double run_seq_wrapper(T *lyst, int NUM, int TLEVEL) {
  quicksort(lyst, NUM);
}

double run_par_wrapper(T *lyst, int NUM, int TLEVEL) {
  parallelQuicksort(lyst, NUM, TLEVEL);
}

double run_libc_wrapper(T *lyst, int NUM, int TLEVEL) {
  qsort(lyst, NUM, sizeof(T), compare_doubles);
}

/*
Main method:
-generate random list
-time sequential quicksort
-time parallel quicksort
-time standard qsort
*/
int main(int argc, char *argv[])
{
  struct timer1 timer;

  srand(47);            // fixed seed

  int NUM = DNUM;
  int TLEVEL = THREAD_LEVEL;

  if (argc > 1)                //user specified list size.
  {
    NUM = atoi(argv[1]);
  }

  if (argc > 2)                //user specified thread level.
  {
    TLEVEL = atoi(argv[2]);
  }


  //Want to compare sorting on the same list,
  //so backup.
  T *lystbck = (T *) malloc(NUM * sizeof(T));
  T *lyst = (T *) malloc(NUM * sizeof(T));

  //Populate random original/backup list.
  for (int i = 0; i < NUM; i++) {
    lystbck[i] = 1.0 * rand() / RAND_MAX;
  }

#ifdef PROFILING
    memcpy(lyst, lystbck, NUM * sizeof(T));
    run_par_wrapper(lyst, NUM, TLEVEL);
    exit(0);
#endif


#ifdef ALL_IMPL

  //Init the array that specify the order in which we will
  //call the sorting functions
  int method[3*REP];
  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < REP; j++) {
      method[i*REP+j] = i;
    }
  }
  std::random_shuffle(method, method+3*REP);

  double (*sorting_function[3])(T*, int, int) =
   {run_seq_wrapper, run_par_wrapper, run_libc_wrapper};

  //Holds the average time taken by each method
  double avg[3] = {0,0,0};

  for (size_t i = 0; i < 3*REP; i++) {
    memcpy(lyst, lystbck, NUM * sizeof(T));

    timer.start();
    sorting_function[method[i]](lyst, NUM, TLEVEL);
    timer.end();

    if (!isSorted(lyst, NUM)) {
      printf("Oops, lyst did not get sorted by method %d.\n", method[i]+1);
      exit(1);
    }

    avg[method[i]] += timer.get_elapsed_time();
  }

  for (size_t i = 0; i < 3; i++) {
    avg[i] /= REP;
  }

  printf("%lf %lf %lf\n", avg[0], avg[1], avg[2]);

#else

  double avg = 0;
  for (size_t i = 0; i < REP; i++) {
    memcpy(lyst, lystbck, NUM * sizeof(T));

    timer.start();
    parallelQuicksort(lyst, NUM, TLEVEL);
    timer.end();

    avg += timer.get_elapsed_time();
  }
  avg /= REP;
  printf("%lf\n", avg);

#endif

  free(lyst);
  free(lystbck);
  pthread_exit(NULL);
}

void quicksort(T lyst[], int size)
{
  quicksortHelper(lyst, 0, size - 1);
}

void quicksortHelper(T lyst[], int lo, int hi)
{
  if (lo >= hi)
    return;
  int b = partition(lyst, lo, hi);
  quicksortHelper(lyst, lo, b - 1);
  quicksortHelper(lyst, b + 1, hi);
}

void swap(T lyst[], int i, int j)
{
  T temp = lyst[i];
  lyst[i] = lyst[j];
  lyst[j] = temp;
}

int partition(T lyst[], int lo, int hi)
{
  int b = lo;
  int r = (int) (lo + (hi - lo) * (1.0 * rand() / RAND_MAX));
  T pivot = lyst[r];
  swap(lyst, r, hi);
  for (int i = lo; i < hi; i++) {
    if (lyst[i] < pivot) {
      swap(lyst, i, b);
      b++;
    }
  }
  swap(lyst, hi, b);
  return b;
}


/*
parallel quicksort top level:
instantiate parallelQuicksortHelper thread, and that's
basically it.
*/
void parallelQuicksort(T lyst[], int size, int tlevel)
{
  int rc;
  void *status;

  //Want joinable threads (usually default).
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  //pthread function can take only one argument, so struct.
  struct thread_data td;
  td.lyst = lyst;
  td.low = 0;
  td.high = size - 1;
  td.level = tlevel;

  //The top-level thread.
  pthread_t theThread;
  rc = pthread_create(&theThread, &attr, parallelQuicksortHelper,
                      (void *) &td);
  if (rc) {
    printf("ERROR; return code from pthread_create() is %d\n", rc);
    exit(-1);
  }
  //Now join the thread (wait, as joining blocks) and quit.
  pthread_attr_destroy(&attr);
  rc = pthread_join(theThread, &status);
  if (rc) {
    printf("ERROR; return code from pthread_join() is %d\n", rc);
    exit(-1);
  }
  //printf("Main: completed join with top thread having a status of %ld\n",
  //              (long)status);

}

/*
parallelQuicksortHelper
-if the level is still > 0, then partition and make
parallelQuicksortHelper threads to solve the left and
right-hand sides, then quit. Otherwise, call sequential.
*/
void *parallelQuicksortHelper(void *threadarg)
{
  int mid, t, rc;
  void *status;

  struct thread_data *my_data;
  my_data = (struct thread_data *) threadarg;

  //fyi:
  //printf("Thread responsible for [%d, %d], level %d.\n",
  //              my_data->low, my_data->high, my_data->level);

  if (my_data->level <= 0 || my_data->low == my_data->high+1) {
    //We have plenty of threads, finish with sequential.
    quicksortHelper(my_data->lyst, my_data->low, my_data->high);
    pthread_exit(NULL);
  }
  //Want joinable threads (usually default).
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  //Now we partition our part of the lyst.
  mid = partition(my_data->lyst, my_data->low, my_data->high);

  //At this point, we will create threads for the
  //left and right sides.  Must create their data args.
  struct thread_data thread_data_array[2];

  for (t = 0; t < 2; t++) {
    thread_data_array[t].lyst = my_data->lyst;
    thread_data_array[t].level = my_data->level - 1;
  }
  thread_data_array[0].low = my_data->low;
  thread_data_array[0].high = mid - 1;
  thread_data_array[1].low = mid + 1;
  thread_data_array[1].high = my_data->high;

  //Now, instantiate the threads.
  //In quicksort of course, due to the transitive property,
  //no elements in the left and right sides of mid will have
  //to be compared again.
  pthread_t threads[2];
  for (t = 0; t < 2; t++) {
    rc = pthread_create(&threads[t], &attr, parallelQuicksortHelper,
                        (void *) &thread_data_array[t]);
    if (rc) {
      printf("ERROR; return code from pthread_create() is %d\n", rc);
      exit(-1);
    }
  }

  pthread_attr_destroy(&attr);
  //Now, join the left and right sides to finish.
  for (t = 0; t < 2; t++) {
    rc = pthread_join(threads[t], &status);
    if (rc) {
      printf("ERROR; return code from pthread_join() is %d\n", rc);
      exit(-1);
    }
  }

  pthread_exit(NULL);
}

//check if the elements of lyst are in non-decreasing order.
//one is success.
int isSorted(T lyst[], int size)
{
  for (int i = 1; i < size; i++) {
    if (lyst[i] < lyst[i - 1]) {
      printf("at loc %d, %e < %e \n", i, lyst[i], lyst[i - 1]);
      return 0;
    }
  }
  return 1;
}

//for the built-in qsort comparator
//from http://www.gnu.org/software/libc/manual/html_node/Comparison-Functions.html#Comparison-Functions
int compare_doubles(const void *a, const void *b)
{
  const T *da = (const T *) a;
  const T *db = (const T *) b;

  return (*da > *db) - (*da < *db);
}
