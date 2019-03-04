#ifndef H_CBICA_UTIL_GTS
#define H_CBICA_UTIL_GTS

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