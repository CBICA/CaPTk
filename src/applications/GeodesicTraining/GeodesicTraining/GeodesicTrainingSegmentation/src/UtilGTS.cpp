#include "UtilGTS.h"

std::string GeodesicTrainingSegmentation::UtilGTS::currentDateTime()
{
	// Get current time
	std::chrono::time_point<std::chrono::system_clock> time_now = std::chrono::system_clock::now();

	// Convert to time_t for <ctime>
	std::time_t time_now_t = std::chrono::system_clock::to_time_t(time_now);

	// Format to datetime
	std::tm now_tm = *std::localtime(&time_now_t);
	char buf[512];
	std::strftime(buf, 512, "%Y-%m-%d %H.%M.%S", &now_tm);
	std::string dateTime(buf);

	return dateTime;
}