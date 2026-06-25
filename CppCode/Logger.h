#ifndef LOGGER_H_INCLUDED
#define LOGGER_H_INCLUDED

#include <iostream>
#include <iomanip>
#include <vector>

// @brief Wrapper class for writing data (particularly vectors) to a text file
// @app Boolean to determine if data should be appended to file
// @fname String indicating location and name of file (full path)
template <typename T>
struct Logger
{
	Logger() = default;
	Logger(const std::string& fname, const bool& isAppendable) 
	{
	    if (isAppendable) 
	        Stream.open(fname.c_str(), std::ios::out | std::ios::app);
	    else 
	        Stream.open(fname.c_str(), std::ios::out);
	    Stream << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;
	}
	std::ofstream Stream;
};

//template <typename T>
//struct Loggers
//{
//	Loggers(const std::vector<std::string>& fname, const bool& isAppendable) 
//	{
//	    Logs.resize(fname.size());
//	    for (size_t i = 0; i < fname.size(); ++i) {
//		Logs[i] = Logger<T>(fname[i], isAppendable);
//		Logs[i] << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;
//	    }
//	}
//	std::vector<Logger<T>> Logs;
//};


#endif
