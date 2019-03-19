#include "SvmSuiteUtil.h"

SvmSuiteUtil::Timer::Timer() {
	Reset();
}

void SvmSuiteUtil::Timer::Reset() {
	m_timestamp = std::chrono::high_resolution_clock::now();
}

float SvmSuiteUtil::Timer::Diff() {
	std::chrono::duration<float> fs = std::chrono::high_resolution_clock::now() - m_timestamp;
	return fs.count();
}