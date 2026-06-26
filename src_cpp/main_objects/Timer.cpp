#include "Timer.h"

#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>

// PRIVATE METHODS

std::string Timer::formatElapsed(double seconds) const
{
    long long total = static_cast<long long>(seconds);

    long long hours   = total / 3600;
    long long minutes = (total % 3600) / 60;
    long long secs    = total % 60;

    std::ostringstream out;
    out << std::setfill('0')
        << std::setw(2) << hours   << ":"
        << std::setw(2) << minutes << ":"
        << std::setw(2) << secs;

    return out.str();
}

// PUBLIC METHODS

void Timer::start()
{
    start_steady = std::chrono::steady_clock::now();
    start_wall   = std::chrono::system_clock::now();
}

double Timer::elapsedSeconds() const
{
    auto now = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(now - start_steady).count();
}

void Timer::printInterval(const std::string& label) const
{
    std::cout << label
              << formatElapsed(elapsedSeconds());
}

void Timer::stop() const
{
    auto end_steady = std::chrono::steady_clock::now();
    auto end_wall   = std::chrono::system_clock::now();

    double elapsed =
        std::chrono::duration<double>(end_steady - start_steady).count();

    std::time_t start_time =
        std::chrono::system_clock::to_time_t(start_wall);

    std::time_t end_time =
        std::chrono::system_clock::to_time_t(end_wall);

    std::cout << "\nStarted : "
              << std::put_time(std::localtime(&start_time),
                               "%Y-%m-%d %H:%M:%S")
              << '\n';

    std::cout << "Finished: "
              << std::put_time(std::localtime(&end_time),
                               "%Y-%m-%d %H:%M:%S")
              << '\n';

    std::cout << "Elapsed : "
            << formatElapsed(elapsed)
            << '\n';

    std::cout << "\nProgram run complete.\n";
}