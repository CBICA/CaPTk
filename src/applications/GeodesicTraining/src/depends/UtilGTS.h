#ifndef H_CBICA_UTIL_GTS
#define H_CBICA_UTIL_GTS

#pragma warning(disable : 4996) //_CRT_SECURE_NO_WARNINGS

#include <ctime>
#include <chrono>
#include <iostream>

namespace GeodesicTrainingSegmentation 
{
	namespace UtilGTS
	{
		/**
		Get current date and time
		@return current date and time as %Y-%m-%d %H.%M.%S
		*/
		std::string currentDateTime();
	}

}

#endif