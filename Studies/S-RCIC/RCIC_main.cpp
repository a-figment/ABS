#include <iostream>
#include <algorithm>
#include <chrono>
#include <thread>
#include "Billiard.h"
#include "Trajectory.h"
#include "OnlineStatistics.h"
#include "Misc.h"
#include "Ensemble.h"
#include "Simulation.h"
#include "DataKeys.h"

// < single character prefix, max number of instantiations >
typedef float_ldub float_;
const int nDigits = std::numeric_limits<float_>::digits10;

int main(int argc, char *argv[])
{
    std::cout << "Arugments received: ";
    for (int i = 0; i < argc; ++i) 
        std::cout << argv[i] << " ";
    std::cout << std::endl;

    // ======================= BOILER PLATE START ==========================
    std::cout << std::setprecision(std::numeric_limits<float_>::digits10) << std::showpoint;
    std::vector<int> KeyIndices;
    for (int i = 2; i < argc; ++i) 
        KeyIndices.push_back(std::atoi(argv[i]));
    DataKeys DK({"p", "b", "d"}, {100, 6, 2}, {"Partitions", "Boundary parameter: [0-2] Irrational, [3-5] Rational", "Pore widths: [0] Finite horizon, [1] Infinite horizon"});
    std::string Key = DK.MakeKey(KeyIndices);
    DK.Print();
    // ======================= BOILER PLATE END ==========================
    
    // ======================= BILLIARD DEFINITION START ==========================
    std::vector<std::unique_ptr<Obstacle<float_>>> Obstacles(1);
    
    // Internal angle definitions
    float_ PHI = boost::math::constants::phi<float_>();
    float_ PI = boost::math::constants::pi<float_>();

    // Contained within [0, pi / 2]
    // std::vector<float_> alpha_irrational = seq((5 * PI) / (2 + 5 * PHI), (5 * PI) / (25 * PHI + 10), 3); // IRRATIONAL
    // std::vector<float_> alpha_rational = seq(PI / (float_)2, PI / (float_)10, 3); // RATIONAL

    // Contained within [0, pi]
    std::vector<float_> alpha_rational = { 2 * PI / ((float_)3), PI / ((float_)2), PI /((float_)3) };
    std::vector<float_> alpha_irrational(alpha_rational.size());
    for (size_t i = 0; i < alpha_irrational.size(); ++i) 
        alpha_irrational[i] = (alpha_rational[i] * (PHI - (3 / ((float_)(5)))));

    // 0-2 IRRATIONAL, 3-6 RATIONAL
    std::vector<float_> alpha = alpha_irrational;
    alpha.insert(alpha.end(), alpha_rational.begin(), alpha_rational.end());

    std::vector<float_> porewidths = { - 1/((float_)2), (float_)1 };

    Polygon<float_> Channel = getParallelChannel<float_>(alpha[KeyIndices[1]], (1 / (float_)2) , porewidths[KeyIndices[2]]);

    
    Channel.Write("Data/" + Key + "BilliardTable.dat");
    Obstacles[0] = std::make_unique<Polygon<float_>>(Channel);
    Billiard<float_> Table(Obstacles);
    // ======================= BILLIARD DEFINITION END ==========================
    
    
    float_ T = (float_)1e5;
    //std::vector<float_> Time = seq<float_>(1/((float_)10), T,  1000);
    //std::vector<float_> Time = getExponentialTimes<float_>(T, 1000);

    Simulation<float_> Sim {
        .Time = {0.0, 10.0, 100.0, 1000.0, 10000.0, 25000.0, 50000.0, 75000.0, 100000.0},
        .L = {0},
        .Symbols = {""},
        .DataDestination = "Data/" + Key,
        .RunTime = ".",
        .MomentStep = (float_)0.5,
        .TotalTime = 100001.0,
        .nInitialConditions = 10000,
        .nCollisions = 0,
        .nMoments = 20};
    Sim.Write();
    
    //WriteSampleTrajectory<float_>(Table, Sim);
    //VelocityCorrelations<float_>(Table, Sim);  
    ContinuousSystemSummary<float_>(Table, Sim);  


    //Simulation<float_> FiniteSim {
    //    .Time = Time,
    //    .L = {(float_)10},
    //    .Symbols = {""},
    //    .DataDestination = "Data/" + Key,
    //    .RunTime = ".",
    //    .MomentStep = 1,
    //    .TotalTime = (float_)Time.back(),
    //    .nInitialConditions = 100000,
    //    .nCollisions = 0,
    //    .nMoments = 10};

    //FiniteSim.Write();
    ////EscapeToAbsorbingBoundary<float_>(Table, FiniteSim);  
    //FirstPassageTimes<float_>(Table, FiniteSim);
    //ForcedEscapeToAbsorbingBoundary<float_>(Table, FiniteSim);  
    //MeanFirstPassage(Table, FiniteSim);


    //Simulation<float_> FiniteSim {
    //    .Time = {0},
    //    .L = {((float_)50)},
    //    .Symbols = {""},
    //    .DataDestination = "Data/" + Key,
    //    .RunTime = ".",
    //    .MomentStep = 0,
    //    .TotalTime = (float_)1e8,
    //    .nInitialConditions = 10000,
    //    .nCollisions = 0,
    //    .nMoments = 0};

    //FiniteSim.Write();
    //FirstPassageTimes(Table, FiniteSim);



    //std::thread t1(ContinuousSystemSummary<float_>, std::ref(Table), std::ref(Sim));
    //std::thread t2(MeanFirstPassage<float_>, std::ref(Table), std::ref(FiniteSim));

    //t1.join();
    //t2.join();

   return 0;
}

