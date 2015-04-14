#define HADGEN_LIB_INTERNAL
#include "hadgen.h"
#include "hadgen_common_blocks.h"
#include <time.h>
#include <signal.h>
#include <pthread.h>
#include <stdio.h>

#if _POSIX_C_SOURCE < 199309L
   #error Current POSIX source version is not supported
#endif

// this file implements a kind of crutch for HADGEN original code
// if generator hangs, this will reset it after timeout

static int timeout = 3;
static pthread_t thread;

void hadgen_set_timeout(int sec) {
   if (sec < 0) sec = -sec;
   if (sec < 2) sec = 2;
   timeout = sec;
}

int hadgen_get_timeout() {
   return timeout;
}

void reset_thread(union sigval sig) {
   pthread_cancel(thread);
}

void init_timer(timer_t timer, struct itimerspec *t_spec) {
   (*t_spec).it_interval.tv_sec = 0;
   (*t_spec).it_interval.tv_nsec = 0;
   (*t_spec).it_value.tv_sec = timeout;
   (*t_spec).it_value.tv_nsec = 0;
   timer_settime(timer, 0, t_spec, 0);
}

void cancel_timer(timer_t timer) {
   struct itimerspec t_spec;
   t_spec.it_interval.tv_sec = 0;
   t_spec.it_interval.tv_nsec = 0;
   t_spec.it_value.tv_sec = 0;
   t_spec.it_value.tv_nsec = 0;
   timer_settime(timer, 0, &t_spec, 0);
}

int hadgen_start() {
   // init thread attributes
   // make it joinable
   pthread_attr_t attr;
   pthread_attr_init(&attr);
   pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
   // status variable
   void *status = 0;
   // create timer
   struct sigevent sevp;
   timer_t timer;
   sevp.sigev_notify = SIGEV_THREAD;
   sevp.sigev_notify_function = reset_thread;
   sevp.sigev_value.sival_int = 0;
   timer_create(CLOCK_REALTIME, &sevp, &timer);
   // timer structures
   struct itimerspec t_spec;
   int first_start = 1;
   do {
      if ((long)status == 1) {
         // Bad parameters on initialization, cannot restart
//         printf("HADGEN: bad parameters, cannot start\n");         
         return 1;
      }
      if (!first_start) {
         // change the random seed
         time_t t = time(0);
         srand(t);
         unsigned long seed = rand();
         hadgen_set_randomseed(seed);
      }
      // start the computation
      pthread_create(&thread, &attr, hadgen_thread, 0);
      // arm the timer
      init_timer(timer, &t_spec);
      // wait for the thread 
      pthread_join(thread, &status);
      // disarm the timer
      cancel_timer(timer);
      first_start = 0;
   } while ((long)status != 101);  // while not success
   pthread_attr_destroy(&attr);
   return 0;
}

void hadgen_cleanup(void *arg) {
   hadgen_terminate();
}

void *hadgen_thread(void *p) {
   if (hadgen_initialize() != 0) pthread_exit((void*)1);
   int ss, ss2;
   pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, &ss);
   pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, &ss); 
   // setup cleanup
   pthread_cleanup_push(hadgen_cleanup, 0);
   // start hadgen
   hadgen_short();
   // exit normally (101 - success code)]
   pthread_cleanup_pop(0);
   pthread_exit((void*)101);
}


