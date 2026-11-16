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
const float_ PHI = boost::math::constants::phi<float_>();
const float_ PI = boost::math::constants::pi<float_>();


// @brief Writes moments of displacement and displacements at selected times
template <typename T>
inline void getMoments_t(Billiard<T> &Table, Simulation<T> &Sim) 
{
    // Involves all initial conditions
    std::vector<Ensemble<T>> Moments(Sim.nMoments, Ensemble<T>(Sim.Time.size()));

    Writer<T> w_DX(Sim.DataDestination + "DX.dat", true);
    std::vector<size_t> Indices = getIndicesOfVectorValue(Sim.TimesToGet, Sim.Time);
    std::cout << "Collecting displacements at times t = " << index(Sim.Time, Indices) << " with indices = " << Indices << std::endl;


    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < Sim.nInitialConditions; ++i)
    {
	if (10 * i % Sim.nInitialConditions == 0 && i != 0)
	    std::cout << "Simulated " << i << " initial conditions." << std::endl;

	Particle<T> P0 = Table.getInitialState(0, -1, Sim.TotalTime);

	Trajectory<T> Traj = Table.getTrajectory(P0, Sim.Time);
        std::vector<T> xDis = Traj.getDisplacements_x();
	w_DX.WriteRowVector(xDis, Indices);

	std::vector<std::vector<T>> DXp = Traj.getDisplacements_x(Sim.MomentStep, Sim.nMoments);
	for (int j = 0; j < Sim.nMoments; ++j)
	    Moments[j].Push(DXp[j]);

    }
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";

    for (size_t i = 0; i < Moments.size(); ++i)
	Moments[i].Write_Append(Sim.DataDestination + "Moments.dat");
}


int main(int argc, char *argv[])
{
    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();

    // DEFINING DATA KEYS
    std::cout << "Arugments received: ";
    for (int i = 0; i < argc; ++i)
        std::cout << argv[i] << " ";
    std::cout << std::endl;
    std::cout << std::setprecision(std::numeric_limits<float_>::digits10) << std::showpoint;
    std::vector<int> KeyIndices;
    for (int i = 2; i < argc; ++i)
        KeyIndices.push_back(std::atoi(argv[i]));

    // DATA KEYS FOR PARALLEL CHANNELS
    DataKeys DK({"p", "a", "r", "t", "n"}, {1000, 1000, 1000, 3, 3}, {"PARTITIONS", "ALPHA", "RATIO", "TOTAL TIME", "NUM. TIME STEPS"});
    std::string Key = DK.MakeKey(KeyIndices);
    DK.Print();

    // A SINGLE OBSTACLE BILLIARD
    std::vector<std::unique_ptr<Obstacle<float_>>> Obstacles(1);

    // RATIONAL CHANNELS DEFINED BY APERTURE ANGLE AND RATIO
    std::vector<float_> alpha_rational = getRationalMultipleOfPi<float_>(7);
    std::vector<float_> alpha(alpha_rational.size() * 2);
    for (size_t i = 0; i < alpha.size(); ++i) {
        if (i < (alpha_rational.size()))
            alpha[i] = alpha_rational[i];
        else
            alpha[i] = (alpha_rational[i - alpha_rational.size()] * (PHI - (3 / ((float_)(5)))));
    }
    float_ dx = (1/((float_)2));
    float_ dy = dx * cos(alpha[KeyIndices[1]] / (float_)2) / sin(alpha[KeyIndices[1]] / (float_)2);

    // Pore width in terms of ratio
    // r[0-30]   -> Finite horizon. { r[0], r[10], r[20], r[30] } = { 1/4, 1/2, 3/4, 1 }
    // r[31-390] -> Ininite horizon. Steps of 20 = increments of 1/2 
    //std::vector<float_> ratio = seq<float_>(1/(float_)5, 10, 393);  
    std::vector<float_> ratio = seq<float_>(1/(float_)5, 10, 981); // [0-80] FH
    std::vector<float_> porewidths(ratio.size(), 0); 
    for (size_t i = 0; i < porewidths.size(); ++i) 
        porewidths[i] = ratio[i] * dy;

    // Create channel
    Polygon<float_> Channel = getParallelChannel2<float_>(dy, (1 / (float_)2), porewidths[KeyIndices[2]]);

    // Summarise parameters
    std::cout << std::setprecision(4) << "(aperture, height, ratio, width) = (" << alpha[KeyIndices[1]] * 180/PI << ", " << dy << ", " << ratio[KeyIndices[2]] << ", " << porewidths[KeyIndices[2]] << ")." << std::endl;


    // NON-PARALLEL CHANNELS
    //std::vector<float_> H = { (float_)2, 3/(float_)2, (float_)1, 3/(float_)4, 1/(float_)2 };
    //float_ dyt = dx;
    //std::vector<float_> dyb_on_dx = {(float_)1 + 1/(float_)100, 
    //    			     (float_)1 + 1/(float_)1000, 
    //    			     (float_)1 + 1/(float_)10000, 
    //    			     (float_)1 + 1/(float_)100000, 
    //    			     (float_)1 + 1/(float_)1000000, 
    //    			     (float_)1}; // = dyb / dx 
    //std::vector<float_> dyb(dyb_on_dx.size());
    //for (size_t i = 0; i < dyb.size(); ++i) 
    //    dyb[i] = dyb_on_dx[i] * dx;
    //Polygon<float_> Channel = getNonParallelChannel<float_>(1/(float_)2, dyb[KeyIndices[2]], dyt, H[KeyIndices[1]] - dyt);

    Obstacles[0] = std::make_unique<Polygon<float_>>(Channel);
    Billiard<float_> Table(Obstacles);


    // PARAMETERS FOR SIMULATION OBJECT
    std::vector<float_> T = { 1e1, 1e2, 1e3 }; // Total time
    std::vector<uint> nSteps = { 101, 1001 };  // Number of time steps

    std::vector<float_> TIME = seq<float_>(0, T[KeyIndices[3]], nSteps[KeyIndices[4]]); // Time
    std::vector<float_> ttg = { TIME[1], TIME[10], TIME[100] };	// Times to inspect an observable
    int M = 1000000;	      // Number of initial conditions
    int nm = 2;		      // Number of moments
    float_ mstep = (float_)1; // Moment step
    
    std::string DataDir = "Data/" + Key;
    //std::string DataDir = "/fred/oz055/pathtohpc/ABS/S-Moments/" + Key;
    
    Simulation<float_> Sim {
        .Time = TIME,
        .L = {0},
        .TimesToGet = ttg,
        .Symbols = {""},
        .DataDestination = DataDir,
        .RunTime = ".",
        .MomentStep = mstep,
        .TotalTime = TIME.back(),
        .nInitialConditions = M,
        .nCollisions = 0, // always reach max time before max coll
        .nMoments = nm};
    Sim.Print();

    getMoments_t<float_>(Table, Sim);

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";
    return 0;
}


    

