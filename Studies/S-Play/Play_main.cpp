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
#include "EscapeData.h"
#include "CoarseGrain.h"

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
    DataKeys DK({"p", "a", "d"}, {1000, 20, 40}, {"PARTITIONS", "ALPHA_R = PI/2", "PORE WIDTH"});
    std::string Key = DK.MakeKey(KeyIndices);
    DK.Print();

    // Consts
    float_ PHI = boost::math::constants::phi<float_>();
    float_ PI = boost::math::constants::pi<float_>();
    // ======================= BOILER PLATE END ==========================
    
    // ======================= BILLIARD DEFINITION START ==========================
    std::vector<std::unique_ptr<Obstacle<float_>>> Obstacles(1);
    std::vector<float_> alpha_rational = getRationalMultipleOfPi<float_>(1); 
    std::vector<float_> alpha(alpha_rational.size() * 2);
    for (size_t i = 0; i < alpha.size(); ++i) {
        if (i < (alpha_rational.size()))
            alpha[i] = alpha_rational[i];
        else
            alpha[i] = (alpha_rational[i - alpha_rational.size()] * (PHI - (3 / ((float_)(5)))));
    }
    float_ dx = (1/((float_)2));
    float_ dy = dx * cos(alpha[KeyIndices[1]] / (float_)2) / sin(alpha[KeyIndices[1]] / (float_)2);
    std::vector<float_> ratio = seq<float_>(1/(float_)4, 10, 40);
    std::vector<float_> porewidths(ratio.size(), 0); //= { - 1/((float_)2), (float_)1 };
    for (size_t i = 0; i < porewidths.size(); ++i) {
        porewidths[i] = ratio[i] * dy;
    }

    Polygon<float_> Channel = getParallelChannel2<float_>(dy, (1 / (float_)2), porewidths[KeyIndices[2]]);
    Obstacles[0] = std::make_unique<Polygon<float_>>(Channel);
    Billiard<float_> Table(Obstacles);
    // ======================= BILLIARD DEFINITION END ==========================

    // =============================== SIMULATION START ===============================    
    std::vector<float_> TIME = seq<float_>(0, 1e6, 1001);
    std::string DataDir = "Data/" + Key;
    int M = 3;
    int N = 10000;
    uint w = 10;


    //Particle<float_> P0 = Table.getInitialState(0, -1, TIME.back());
    //Table.WriteLongTimeSeries(P0, N, DataDir);

    for (int k = 0; k < 1; ++k) {
	Particle<float_> P0 = Table.getInitialState(0, -1, TIME.back());

    	//Trajectory<float_> Traj = Table.getTrajectory(P0, TIME);
    	Trajectory<float_> Traj = Table.getTrajectory(P0, N);
    	std::vector<float_> x = Traj.getLiftedPositions_x();

	CoarseGrain<float_> cg(TIME, x, w);
	cg.WriteTurningPoints(DataDir);
	
	if (cg.isGrainable()) {
	    std::vector<float_> tau = cg.getJumpTime();
	    WriteColVector<float_>(DataDir + "x.dat", x);
	    WriteColVector<float_>(DataDir + "t.dat", Traj.getLiftedFreePaths());
	    WriteVector<float_>(DataDir + "tau.dat", tau);
	}
    	
    }

    //Writer<float_> w_TAU(DataDir + "FP.dat", true);
    //Writer<float_> w_TH(DataDir + "Theta.dat", true);
    //Writer<float_> w_X(DataDir + "EPx.dat", true);
    //Writer<float_> w_Y(DataDir + "EPy.dat", true);
    //Writer<float_> w_THE(DataDir + "EscapeVel.dat", true);
    //Writer<float_> w_ET(DataDir + "EscapeTime.dat", true);

    ////int m = 10;
    ////int n = 10;
    ////std::vector<float_> Bp = seq<float_>(-PI/(float_)2, PI/(float_)2, m); // B^-
    ////std::vector<float_> Bm1 = seq<float_>(PI/(float_)2, PI, m/2); // B^+
    ////std::vector<float_> Bm2 = seq<float_>(-PI, -PI/(float_)2, m/2); // B^+
    ////std::vector<float_> Bm = Bm1;
    ////Bm.insert(Bm.end(), Bm2.begin(), Bm2.end());
    ////std::vector<float_> D = seq<float_>(0, d[KeyIndices[2]], n);  // D

    //std::vector<float_> L = seq(1/(float_)4, 10000 + 1/(float_)2, 10001);

    //for (int i = 0; i < M; ++i) // Loop over M particles
    //{
    //    if (10 * i % M == 0 && i != 0)
    //        std::cout << "Simulated " << i << " initial conditions." << std::endl;

    //    Particle<float_> P0 = Table.getInitialState(0, -1, TIME.back());

    //    Trajectory<float_> Traj = Table.getTrajectory(P0, N);
    //    EscapeData<float_> ED(L, Traj);
    //    ED.setEscapeData(L);

    //    w_Y.WriteRowVector(ED.EscapePosition_y);
    //    w_X.WriteRowVector(ED.EscapePosition_x);
    //    w_THE.WriteRowVector(ED.EscapeDirection);
    //    w_ET.WriteRowVector(ED.EscapeTimes);
    //    
    //    WriteVector<float_>(DataDir + std::to_string(i) + "--xTraj.dat", Traj.getLiftedPositions_x());
    //    WriteVector<float_>(DataDir + std::to_string(i) + "--yTraj.dat", Traj.getLiftedPositions_y());
    //    //Trajectory<float_> Traj = Table.getTrajectory(P0, Sim.Time);
    //    
    //    // Get displacements or other observables from trajectory
    //    //std::vector<float_> tau = Traj.getLiftedFreePaths();
    //    //w_TAU.WriteRowVector(tau);
    //    //w_TH.WriteRowVector(Traj.getVelocities_a());

    //}

    // =============================== SIMULATION END  ===============================    


    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";
    return 0;
}


