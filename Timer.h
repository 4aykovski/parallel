#pragma once
#include <chrono>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <string>

class Timer {
  std::chrono::high_resolution_clock::time_point start;

public:
  Timer();
  double stop();
  std::string get();
  std::string getCurrentTime();
};
