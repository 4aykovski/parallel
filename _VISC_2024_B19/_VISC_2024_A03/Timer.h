#pragma once
#include <ctime>
#include <chrono>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>

class Timer {
	std::chrono::high_resolution_clock::time_point start;
public:
	Timer();
	double stop();
	std::string get();
	std::string getCurrentTime();
};
