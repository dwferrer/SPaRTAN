#include <omp.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <cassert>
#include <stdio.h>

#ifndef _STOPWATCH_CC
#define _STOPWATCH_CC

class StopWatch {
public:
    StopWatch(char *s);
    ~StopWatch();
    void Start(void);
    void Stop(void);
    double Elapsed(void);
    void Clear(void);

    void Print(void) { printf("%20s took %e seconds \n", _name, Elapsed() ); }

    char _name[1024];

    int nprocs;

    int *timeron;

    struct timeval *tuse, *tstart, *timer;
};

StopWatch::StopWatch(char *str) { 

    sprintf(_name,"%s",str);
    nprocs = omp_get_num_procs();
    timeron = false;
    tuse = new timeval[nprocs];
    tstart = new timeval[nprocs];
    timer = new timeval[nprocs];
    timeron = new int[nprocs];

    for(int g=0; g<nprocs; g++) timerclear(&timer[g]);
    for(int g=0; g<nprocs; g++) timeron[g] = 0;
}

StopWatch::~StopWatch() {
    delete[] tuse;
    delete[] tstart;
    delete[] timer;
    delete[] timeron;
}

void StopWatch::Start() {
    int g = omp_get_thread_num();

    if( timeron[g] ) {
        std::cerr << "StopWatch: '" << _name << "' in operation; cannot start!\n";
        exit(1);
    }

    assert( gettimeofday( &(tstart[g]), (struct timezone *)NULL ) == 0 );
    timeron[g] = 1;
}

void StopWatch::Stop(void) {
    int g = omp_get_thread_num();

    if( !timeron[g] ) {
        std::cerr << "StopWatch: '" << _name << "' not in operation; cannot stop!\n";
        exit(1);
    }
    struct timeval dt;
    assert( gettimeofday( &(tuse[g]), (struct timezone *)NULL ) == 0 );
    timersub(&(tuse[g]), &(tstart[g]), &dt);
    timeradd(&dt, &(timer[g]), &(timer[g]));
    timeron[g] = 0;
}

void StopWatch::Clear(void) {
    int sum = 0;
    for(int g=0; g<nprocs; g++) { assert(!timeron[g]); sum += timeron[g]; }

    if(sum!=0) {
        std::cerr << "StopWatch: '" << _name << "' in operation; cannot clear!\n";
        exit(1);
    }
    for(int g=0; g<nprocs; g++) timerclear(&(timer[g]));
}

double StopWatch::Elapsed(void) {
    double sum = 0;
    for(int g=0; g<nprocs; g++) sum += timer[g].tv_sec + 1e-6*timer[g].tv_usec;
    return sum;
}
#endif
