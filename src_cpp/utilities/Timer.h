#ifndef TIMER_H
#define TIMER_H

#include "MPIContext.h"

#include <chrono>
#include <string>

using namespace std;

class Timer {
private:
    const MPIContext& _mpi;    
    chrono::steady_clock::time_point _startSteady;
    chrono::system_clock::time_point _startWall;

    string formatElapsed(double _seconds) const;
    
public:
    explicit Timer(const MPIContext& mpi);

    void start();
    double elapsedSeconds() const;
    void printInterval(const string& label = "") const;
    void printIterationInterval() const;
    void stop() const;
};

#endif