#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <string>

using namespace std;

class Timer {
private:
    chrono::steady_clock::time_point _startSteady;
    chrono::system_clock::time_point _startWall;
    string formatElapsed(double _seconds) const;
    
public:
    void start();
    double elapsedSeconds() const;
    void printInterval(const string& label = "") const;
    void stop() const;
};

#endif