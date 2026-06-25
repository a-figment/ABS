#ifndef ENSEMBLE_INCLUDED
#define ENSEMBLE_INCLUDED

#include "OnlineStatistics.h"
#include "Writer.h"
//#include <vector>
//#include <algorithm>
//#include <iostream>
//#include <fstream>
//#include <string>

template <typename T>
class Ensemble
{
    public:
	Ensemble() = default;
	Ensemble(const size_t &nelems) : nElems(nelems)
	{
    	    OS.resize(nElems);
    	    nUpdates = 0;
    	}

	void Push(const T &x, const size_t &idx)
	{
	    OS[idx].Update(x);
	}

	void Push(const T &x)
	{
	    OS[nUpdates].Update(x);
	    nUpdates++;
	}

	void Push(const std::vector<T> &v)
	{
	    for (size_t i = 0; i < v.size(); ++i)
		OS[i].Update(v[i]);
	    nUpdates++;
	}

	void Push(const std::vector<T> &v, const size_t &range)
	{
	    for (size_t i = 0; i < range; ++i)
		OS[i].Update(v[i]);
	    nUpdates++;
	}

	std::vector<T> PullMeans()
	{
	    std::vector<T> Means(OS.size());
	    for (size_t i = 0; i < Means.size(); ++i)
		Means[i] = OS[i].getMean();
	    return Means;
	}

	void PrintMeans()
	{
	    std::cout << "Means = ";
	    for (size_t i = 0; i < OS.size(); ++i)
		std::cout << OS[i].getMean() << " ";
	    std::cout << std::endl;
	}

	std::vector<T> PullVars() const
	{
	    std::vector<T> Vars(OS.size());
	    for (size_t i = 0; i < Vars.size(); ++i)
		Vars[i] = OS[i].getVariance();
	    return Vars;
	}

	std::vector<long int> PullNumElems() const
	{
	    std::vector<long int> Elems(OS.size());
	    for (size_t i = 0; i < Elems.size(); ++i)
		Elems[i] = OS[i].getNumElems();
	    return Elems;
	}

	// If each average consists of same number of elements, return that number
	long int PullUniqueNumElems() const
	{
	    std::vector<long int> Elems = PullNumElems();
	    if (std::adjacent_find(Elems.begin(), Elems.end(), std::not_equal_to<>()) == Elems.end())
		return Elems[0];
	    return -1;
	}

	void Write(const std::string &Destination)
	{
	    Writer<T> w_Means(Destination, false);
	    // DataDir + "MSD-" + std::to_string(nDigits) + "-digits.dat"
	    //std::ofstream Writer(Destination, std::ios::out);
	    //Writer << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;


	    std::vector<T> means = PullMeans();
	    w_Means.WriteRowVector(means);
	    //for (size_t i = 0; i < means.size(); ++i)
	    //    Writer << means[i] << " ";
	    //Writer << '\n';
	}
	void Write_Append(const std::string &Destination)
	{
	    //Writer<T> w_Means(Destination, true);
	    std::vector<T> means = PullMeans();
	    WriteVector<T>(Destination, means);
	}

	size_t nElems;
	size_t getNumElems() const { return OS.size(); }

    private:
	std::vector<OnlineStatistics<T>> OS;
	size_t nUpdates;
};

#endif // ENSEMBLE_INCLUDED
