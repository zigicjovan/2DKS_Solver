#include "Timer.h"

#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace std;

string Timer::formatElapsed(double _seconds) const {
    long long total = static_cast<long long>(_seconds);
    long long hours = total / 3600;
    long long minutes = (total % 3600) / 60;
    long long secs = total % 60;

    ostringstream out;
    out << setfill('0') << setw(2) << hours << ":" << setw(2) << minutes << ":" << setw(2) << secs;
    return out.str();
}

void Timer::start() {
    _startSteady = chrono::steady_clock::now();
    _startWall = chrono::system_clock::now();
}

double Timer::elapsedSeconds() const {
    auto now = chrono::steady_clock::now();
    return chrono::duration<double>(now - _startSteady).count();
}

void Timer::printInterval(const string& label) const {
    cout << label << formatElapsed(elapsedSeconds());
}

void Timer::stop() const {
    auto endSteady = chrono::steady_clock::now();
    auto endWall   = chrono::system_clock::now();
    double elapsed = chrono::duration<double>(endSteady - _startSteady).count();

    time_t startTime = chrono::system_clock::to_time_t(_startWall);
    time_t endTime = chrono::system_clock::to_time_t(endWall);
    cout << "\nStarted : " << put_time(localtime(&startTime), "%Y-%m-%d %H:%M:%S") << '\n';
    cout << "Finished: " << put_time(localtime(&endTime), "%Y-%m-%d %H:%M:%S") << '\n';
    cout << "Elapsed : " << formatElapsed(elapsed) << '\n';
    cout << "\nProgram run complete.\n";
}