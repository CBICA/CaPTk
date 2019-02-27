#ifndef H_CBICA_SVM_SUITE_UTIL
#define H_CBICA_SVM_SUITE_UTIL

#include <chrono>
#include <ctime>

namespace SvmSuiteUtil
{
	class Timer {
	public:
		/**
		Timer Constructor
		Also instantiates interval start.
		Use Reset() to reset it again at another time.
		*/
		Timer();

		/**
		Resets the start of the interval
		*/
		void Reset();

		/**
		Calculate the difference between reset (or constructor called) and now
		@return the difference in float seconds
		*/
		float Diff();
	private:
		std::chrono::high_resolution_clock::time_point m_timestamp;
	};

}

#endif