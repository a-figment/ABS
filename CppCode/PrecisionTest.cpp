#include <iostream>
#include <algorithm>
#include <chrono>
#include "Billiard.h"
#include "Trajectory.h"
#include "OnlineStatistics.h"
#include "Misc.h"
#include "Ensemble.h"

typedef float_ldub float_;

struct SimulationDetails
{

    SimulationDetails(const int &M, const float_ &TotalTime)
    {
    }

    void InitialisationSummary() // Summary of initial data
    {
        std::cout << nInitialConditions << " initial conditions." << std::endl;
        std::cout << TotalTime << " time units." << std::endl;
    }

    void Summary() //
    {
        std::cout << nInitialConditions << " initial conditions." << std::endl;
        std::cout << TotalTime << " time units." << std::endl;
        std::cout << runtime << " time units." << std::endl;
    }

    int nInitialConditions;
    float_ TotalTime;
    std::string runtime;
};

const int nDigits = std::numeric_limits<float_>::digits10;

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cout << "Not enough command line arguments." << std::endl;
        std::cout << "Require:    Job ID" << std::endl;
        exit(EXIT_FAILURE);
    }
    const int JobID = std::atoi(argv[1]);
    std::cout << "Running job " << argv[1] << std::endl;
    std::string DataDir = "../Studies/S-PrecisionTest/Data/" + std::to_string(std::atoi(argv[1])) + "-";
    std::cout << "Saving data to: " << DataDir << std::endl;
    std::cout << "No. precision digits: " << std::numeric_limits<float_>::digits10 << std::endl
              << std::endl;
    std::cout << std::setprecision(std::numeric_limits<float_>::digits10) << std::showpoint;

    // ======= BILLIARD DEFINITION START =======
    std::vector<std::unique_ptr<Obstacle<float_>>> Obstacles(1);
    std::vector<float_> s = seq(6 / static_cast<float_>(7), 7 / static_cast<float_>(6), 4);
    reverse(s.begin(), s.end());
    std::vector<float_> alpha = multiply<float_>(s, boost::math::constants::pi<float_>() / boost::math::constants::phi<float_>());

    // int jid = JobID % s.size();
    int jid = JobID;
    Polygon<float_> Channel = getParallelChannel<float_>(alpha[jid], 1 / static_cast<float_>(2), 0);
    Channel.Write(DataDir + "Channel.dat");
    Obstacles[0] = std::make_unique<Polygon<float_>>(Channel);
    Billiard<float_> Table(Obstacles);
    // ======= BILLIARD DEFINITION END =======

    // ===============
    float_ T = 1e4; // Total time
    int M = 10000;  // Initial conditions

    std::vector<float_> Time = getExponentialTimes<float_>(T, 1000);

    Ensemble<float_> MSD(Time.size());

    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < M; ++i)
    {
        if (10 * i % M == 0 && i != 0)
            std::cout << "Simulated " << i << " initial conditions." << std::endl;

        Particle<float_> IC = Table.getInitialState(0, -1, T);
        IC.WriteState(DataDir + "InitialConditions.dat");
        Particle<float_> Pt = IC;

        Trajectory<float_> ContTraj = Table.getContinuousTrajectory(Pt, Time);

        if (Table.Failed)
        {
            Table.ErrorReport(Pt, i);
            i--; // Want M initial conditions
        }

        std::vector<Point<float_>> xp = ContTraj.getPositions();

        std::vector<float_> dx2 = ContTraj.getDisplacements_x2();

        MSD.Push(dx2);
    }

    MSD.Write(DataDir + "MSD-" + std::to_string(nDigits) + "-digits.dat");

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    std::cout << "\033[1;32mSimulation run time: " << std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count() << " hours = " << std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count() << " minutes = " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << " milliseconds. \n\033[0m" << std::endl;

    return 0;
}
