#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <string>

class Timer
{
private:
    std::chrono::steady_clock::time_point start_steady;
    std::chrono::system_clock::time_point start_wall;
    std::string formatElapsed(double seconds) const;
    
public:
    void start();
    double elapsedSeconds() const;
    void printInterval(const std::string& label = "") const;
    void stop() const;
};

#endif