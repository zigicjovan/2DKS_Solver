#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <string>

class Timer {
private:
    std::chrono::steady_clock::time_point _startSteady;
    std::chrono::system_clock::time_point _startWall;
    std::string formatElapsed(double _seconds) const;
    
public:
    void start();
    double elapsedSeconds() const;
    void printInterval(const std::string& label = "") const;
    void stop() const;
};

#endif