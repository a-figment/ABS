#include <iostream>
#include <algorithm>
#include <chrono>
#include <thread>
#include "Billiard.h"
#include "Trajectory.h"
#include "OnlineStatistics.h"
#include "Misc.h"
#include "Helper.h"
#include "Ensemble.h"
#include "Simulation.h"
#include "DataKeys.h"
#include "Writer.h"

// THIS TEMPLATE INTENDS TO RESET THE STANDARD FOR 'STUDIES'. THE GOAL IS THE AVOID
// USING THE MISC.H FUNCTIONS TO RUN BILLIARD SIMULATIONS UNLESS IT IS COMPLETELY GENERIC.
//
// EXAMPLE USAGE OF MISC.H FUNCTIONS: 
//		1. COLLECTING MOMENTS OF DISPLACEMENT
//		2. CORRELATION FUNCTIONS (POSITION, VELOCITY, INTEGRATED VELOCITY ETC.)
//		3. 



// < single character prefix, max number of instantiations >
typedef float_ldub float_;
const int nDigits = std::numeric_limits<float_>::digits10;


int main(int argc, char *argv[])
{
    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();


    // ======================================= BOILER PLATE START =======================================    
    std::cout << "Arugments received: ";
    for (int i = 0; i < argc; ++i) 
        std::cout << argv[i] << " ";
    std::cout << std::endl;

    std::cout << std::setprecision(std::numeric_limits<float_>::digits10) << std::showpoint;
    std::vector<int> KeyIndices;
    for (int i = 2; i < argc; ++i) 
        KeyIndices.push_back(std::atoi(argv[i]));
    DataKeys DK({"p", "a", "r"}, {1000, 1000, 1000}, {"PARTITIONS", "ALPHA_R = PI/4, PI/3, ..., 3PI/4", "RATIO"});
    std::string Key = DK.MakeKey(KeyIndices);
    DK.Print();

    // Consts
    float_ PHI = boost::math::constants::phi<float_>();
    float_ PI = boost::math::constants::pi<float_>();
    // ======================= BOILER PLATE END ==========================
    
    // ======================= BILLIARD DEFINITION START ==========================
    std::vector<float_> alpha_rational = getRationalMultipleOfPi<float_>(5); 
    std::vector<float_> alpha(alpha_rational.size() * 2);
    for (size_t i = 0; i < alpha.size(); ++i) {
        if (i < (alpha_rational.size()))
            alpha[i] = alpha_rational[i];
        else
            alpha[i] = (alpha_rational[i - alpha_rational.size()] * (PHI - (3 / ((float_)(5)))));
    }
    float_ dx = (1/((float_)2));
    float_ dy = dx * cot(alpha[KeyIndices[1]] / (float_)2);
    // ==========================================================

    std::vector<float_> ratio = seq<float_>(1/(float_)4, 10, 40);
    std::vector<float_> porewidths(ratio.size(), 0); //= { - 1/((float_)2), (float_)1 };
    for (size_t i = 0; i < porewidths.size(); ++i) {
        porewidths[i] = ratio[i] * dy;
    }
    Polygon<float_> Channel = getParallelChannel2<float_>(dy, (1 / (float_)2), porewidths[KeyIndices[2]]);
    Obstacles[0] = std::make_unique<Polygon<float_>>(Channel);
    Billiard<float_> Table(Obstacles);
    // ==================== PARAMETER SUMMARY END   ========================
    

    // =============================== SIMULATION START ===============================    
    std::vector<float_> TIME = seq<float_>(0, 1e8+1, 2001);

    // CONSTANTS
    std::string DataDestination = "/fred/oz055/jorchard_data/ABS/SAD/";
    int M = 200; // Num. particles
    float_ MomentStep = 1/(float_)2;
    size_t nMoments = 20;
    //

    // OBSERVABLES
    std::vector<Ensemble<float_>> Moments(nMoments, Ensemble<float_>(TIME.size()));
    //

    // MAIN LOOP
    for (int i = 0; i < M; ++i) 
    {
        if (10 * i % M == 0 && i != 0)
            std::cout << "Simulated " << i << " initial conditions." << std::endl;

        Particle<float_> P0 = Table.getInitialState(0, -1, TIME.back());
        Trajectory<float_> Traj = Table.getTrajectory(P0, TIME);

        // Extract trajectory properties
        std::vector<std::vector<float_>> DXp = Traj.getDisplacements_x(MomentStep, nMoments);
        for (int j = 0; j < nMoments; ++j)
            Moments[j].Push(DXp[j]);
        // Get displacements or other observables from trajectory
        std::vector<float_> xDis = Traj.getDisplacements_x();
    }


    for (size_t i = 0; i < Moments.size(); ++i)
        Moments[i].Write(DataDestination + "Moment_" + std::to_string(i) + ".dat");
    Sim.Print(RunTime);

    // =============================== SIMULATION END  ===============================    



    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";
    return 0;
}


    //std::vector<float_> TIME = seq<float_>(1, 1e7+1, 10001);
    //std::vector<float_> ttg = TIME;
    //ttg.erase(ttg.begin());
    //Simulation<float_> Sim {
    //    .Time = TIME,
    //    .L = {0},
    //    .TimesToGet = ttg,
    //    .Symbols = {""},
    //    .DataDestination = "Data/" + Key,
    //    .RunTime = ".",
    //    .MomentStep = 0,
    //    .TotalTime = TIME.back(),
    //    .nInitialConditions = 2000,
    //    .nCollisions = 1, // always reach max time before max coll
    //    .nMoments = 0};
    //Sim.Write(); // Write the simulation object?
    

