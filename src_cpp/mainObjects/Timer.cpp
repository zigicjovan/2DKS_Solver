#include "Timer.h"

#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>

std::string Timer::formatElapsed(double _seconds) const {
    long long total = static_cast<long long>(_seconds);
    long long hours = total / 3600;
    long long minutes = (total % 3600) / 60;
    long long secs = total % 60;

    std::ostringstream out;
    out << std::setfill('0') << std::setw(2) << hours << ":" << std::setw(2) << minutes << ":" << std::setw(2) << secs;
    return out.str();
}

void Timer::start() {
    _startSteady = std::chrono::steady_clock::now();
    _startWall = std::chrono::system_clock::now();
}

double Timer::elapsedSeconds() const {
    auto now = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(now - _startSteady).count();
}

void Timer::printInterval(const std::string& label) const {
    std::cout << label << formatElapsed(elapsedSeconds());
}

void Timer::stop() const {
    auto endSteady = std::chrono::steady_clock::now();
    auto endWall   = std::chrono::system_clock::now();
    double elapsed = std::chrono::duration<double>(endSteady - _startSteady).count();

    std::time_t startTime = std::chrono::system_clock::to_time_t(_startWall);
    std::time_t endTime = std::chrono::system_clock::to_time_t(endWall);
    std::cout << "\nStarted : " << std::put_time(std::localtime(&startTime), "%Y-%m-%d %H:%M:%S") << '\n';
    std::cout << "Finished: " << std::put_time(std::localtime(&endTime), "%Y-%m-%d %H:%M:%S") << '\n';
    std::cout << "Elapsed : " << formatElapsed(elapsed) << '\n';
    std::cout << "\nProgram run complete.\n";
}