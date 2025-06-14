#define _CRT_SECURE_NO_WARNINGS 
#include "Timer.h"

Timer::Timer() {
	this->start = std::chrono::high_resolution_clock::now();
}
double Timer::stop() {
	std::chrono::high_resolution_clock::time_point end =
		std::chrono::high_resolution_clock::now();
	std::chrono::duration <double> time_span = std::chrono::duration_cast
		<std::chrono::duration<double>>(end - start);
	return time_span.count();
}
std::string Timer::getCurrentTime() {
	std::stringstream ss;
	std::chrono::system_clock::time_point today = std::chrono::system_clock::now();
	time_t mytime = std::chrono::system_clock::to_time_t(today);
//	ss << std::put_time(localtime(&mytime), "%T");
	ss << std::put_time(localtime(&mytime), "%R");
	return ss.str();
}
std::string Timer::get() {
	double end = this->stop();
	size_t h = end / 3600;
	end -= h * 3600;
	size_t m = end / 60;
	end -= m * 60;
	size_t s = std::round(end);
	std::stringstream ss;
	ss << std::setw(2) << std::setfill('0') << h << ':' <<
		std::setw(2) << std::setfill('0') << m << ':' << 
		std::setw(2) << std::setfill('0') << s << "  " << getCurrentTime();
		//std::setw(2) << std::setfill('0') << s << "  -  " << getCurrentTime();
	return ss.str();
}