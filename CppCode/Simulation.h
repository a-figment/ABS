#ifndef SIMULATION_H_INCLUDED
#define SIMULATION_H_INCLUDED

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

template <typename T>
struct Simulation
{
    void Print();
    void Print(const std::string &runtime);
    void Write();

    std::vector<T> Time;
    std::vector<T> L; // For finite domains
    std::vector<T> TimesToGet;

    std::vector<std::string> Symbols;
    std::string DataDestination;
    std::string RunTime;

    T MomentStep;
    T TotalTime;

    int nInitialConditions;
    int nCollisions;
    int nMoments;
};

    template <typename T>
void Simulation<T>::Print()
{
    std::cout << "=========== Simulation Info. ===========\n";
    std::cout << "M = " << nInitialConditions << "\n"
	<< "N = " << nCollisions << "\n"
	<< "Data destination = " << DataDestination << "\n"
	<< "Number of moments = " << nMoments << "\n"
	<< "Moment step = " << std::setprecision(3) << MomentStep << "\n"
	<< "Time = " << std::setprecision(3) << std::fixed << Time[0] << ", " << Time[1] << ", " << Time[2] << ", ...," << Time[Time.size()-2] << ", " << Time.back() << std::endl;

    std::cout << "Times to get = ";
    for (size_t i = 0; i < TimesToGet.size(); ++i)
	if (i < TimesToGet.size()-1)
	    std::cout << TimesToGet[i] << ", ";
	else
	    std::cout << TimesToGet[i] << "." << std::endl;

    std::cout << "========================================" << std::endl << std::endl;
}

//template <typename T>
//void Simulation<T>::Print()
//{
//    std::cout << "Num. Initial conditions = " << nInitialConditions << "\n"
//              << "Num. Collisions = " << nCollisions << "\n"
//              << "Total time units = " << TotalTime << "\n"
//              << "No. precision digits = " << std::numeric_limits<T>::digits10 << std::endl;
//}

template <typename T>
void Simulation<T>::Write()
{
    std::ofstream Writer(DataDestination + "Simulation.txt", std::ios::out);
    Writer << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;
    if (Writer.is_open())
    {

        Writer << "# Continuous times " << std::endl
               << Time << std::endl;

        Writer << "# Number of initial conditions " << std::endl
               << nInitialConditions << std::endl
               << std::endl;

        Writer << "# Number of collisions " << std::endl
               << nCollisions << std::endl
               << std::endl;

        Writer << "# Total time units " << std::endl
               << TotalTime << std::endl
               << std::endl;

        Writer << "# Lengths " << std::endl;
        Writer << L << std::endl;

        Writer << "# Symbols " << std::endl;
        for (auto x : Symbols)
            Writer << x << " ";
        Writer << std::endl;

        Writer << "# Data destination " << std::endl
               << DataDestination << std::endl
               << std::endl;

        Writer << "# Run time " << std::endl
               << RunTime << std::endl;

        Writer << "# No. precision digits " << std::endl
               << std::numeric_limits<T>::digits10 << std::endl;
    }
    else
    {
        std::cout << "Issue opening file: " << DataDestination + "Simulation.txt" << std::endl;
    }
    Writer.close();
}


#endif
