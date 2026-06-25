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
    DataKeys DK({"p", "a", "r"}, {1000, 10, 100}, {"PARTITIONS", "ALPHA_R = PI/4, PI/3, ..., 3PI/4", "RATIO"});
    std::string Key = DK.MakeKey(KeyIndices);
    DK.Print();

    // Consts
    float_ PHI = boost::math::constants::phi<float_>();
    float_ PI = boost::math::constants::pi<float_>();
    // ======================= BOILER PLATE END ==========================
    
    // ======================= BILLIARD DEFINITION START ==========================
    std::vector<std::unique_ptr<Obstacle<float_>>> Obstacles(1);
    //
    //
    //T epsilon_L = 3/(T)10;
    //T epsilon_U = 4/(T)10;
//    Polygon<float_> assymetricChannel({0.0, 0.3, 0.6, 1.0, 1.0, 0.7, 0.4, 0.0},
//				      {0.0, 0.3, 0.4, 0.0, 0.5, 0.8, 0.9, 0.5},
    Polygon<float_> assymetricChannel({-0.5, -0.2, 0.1, 0.5, 0.5, 0.2, -0.1, -0.5},
				      {0.0 ,  0.3, 0.4, 0.0, 0.5, 0.8, 0.9, 0.5},
                       {3, 7}, false);
    Obstacles[0] = std::make_unique<Polygon<float_>>(assymetricChannel);
    Billiard<float_> Table(Obstacles);
    std::cout << assymetricChannel << std::endl;
    // ======================= BILLIARD DEFINITION END ==========================


    // ==================== PARAMETER SUMMARY START ========================
    //
    // ==================== PARAMETER SUMMARY END   ========================
    

    std::vector<float_> TIME = seq<float_>(1, 1e4, 11);
    //std::vector<float_> ttg = {18000., 36000.};
    std::vector<float_> ttg = {0., 4000., 8000.};
    ttg.erase(ttg.begin());
    Simulation<float_> Sim {
        .Time = TIME,
        .L = {0},
        .TimesToGet = ttg,
        .Symbols = {""},
        .DataDestination = "Data/" + Key,
        .RunTime = ".",
        .MomentStep = 0,
        .TotalTime = TIME.back(),
        .nInitialConditions = 1000000,
        .nCollisions = 1, // always reach max time before max coll
        .nMoments = 0};

    // =============================== SIMULATION START ===============================    
    WriteTrajectory_All<float_>(Table, Sim);
    int M = 100;

    Writer<float_> w_DX(Sim.DataDestination + "DX.dat", true);
    std::vector<size_t> Indices = getIndicesOfVectorValue(Sim.TimesToGet, Sim.Time);
    std::cout << "Collecting displacements at times t = " << index(Sim.Time, Indices) << " with indices = " << Indices << std::endl;
    for (int i = 0; i < Sim.nInitialConditions; ++i) // Loop over M particles
    {
        if (10 * i % Sim.nInitialConditions == 0 && i != 0)
            std::cout << "Simulated " << i << " initial conditions." << std::endl;

        Particle<float_> P0 = Table.getInitialState(0, -1, Sim.Time.back());
	//std::cout << P0 << std::endl;

        //Trajectory<T> Traj = Table.getTrajectory(P0, Sim.nCollisions);
        Trajectory<float_> Traj = Table.getTrajectory(P0, Sim.Time);
	
        // Get displacements or other observables from trajectory
        std::vector<float_> xDis = Traj.getDisplacements_x();
	w_DX.WriteRowVector(xDis, Indices);
    }



    // =============================== SIMULATION END  ===============================    


    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";
    return 0;
}


