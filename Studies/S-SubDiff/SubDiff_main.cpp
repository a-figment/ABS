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
//typedef float_50 float_;
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
    DataKeys DK({"p", "h", "y"}, {100, 4, 6}, {"Partitions", "Height", "dy_b/dx"});
    std::string Key = DK.MakeKey(KeyIndices);
    DK.Print();
    // ======================= BOILER PLATE END ==========================
    
    // ======================= BILLIARD DEFINITION START ==========================
    std::vector<std::unique_ptr<Obstacle<float_>>> Obstacles(1);
    
    // Internal angle definitions
    //float_ PHI = boost::math::constants::phi<float_>();
    //float_ PI = boost::math::constants::pi<float_>();

    //// Contained within [0, pi / 2]
    //// std::vector<float_> alpha_irrational = seq((5 * PI) / (2 + 5 * PHI), (5 * PI) / (25 * PHI + 10), 3); // IRRATIONAL
    //// std::vector<float_> alpha_rational = seq(PI / (float_)2, PI / (float_)10, 3); // RATIONAL

    //// Contained within [0, pi]
    //std::vector<float_> alpha_rational = { 2 * PI / ((float_)3), PI / ((float_)2), PI /((float_)3) };
    //std::vector<float_> alpha_irrational(alpha_rational.size());
    //for (size_t i = 0; i < alpha_irrational.size(); ++i) 
    //    alpha_irrational[i] = (alpha_rational[i] * (PHI - (3 / ((float_)(5)))));

    //// 0-2 IRRATIONAL, 3-6 RATIONAL
    //std::vector<float_> alpha = alpha_irrational;
    //alpha.insert(alpha.end(), alpha_rational.begin(), alpha_rational.end());

    //std::vector<float_> porewidths = { - 1/((float_)2), (float_)1 };

    float_ dx = 1/(float_)2;
    std::vector<float_> H = { (float_)2, 3/(float_)2, (float_)1, 3/(float_)4 };
    float_ dyt = dx;
    std::vector<float_> dyb_on_dx = {(float_)1 + 1/(float_)100, 
				     (float_)1 + 1/(float_)1000, 
				     (float_)1 + 1/(float_)10000, 
				     (float_)1 + 1/(float_)100000, 
				     (float_)1 + 1/(float_)1000000, 
				     (float_)1}; // = dyb / dx 

    std::vector<float_> dyb(dyb_on_dx.size());
    for (size_t i = 0; i < dyb.size(); ++i) {
	dyb[i] = dyb_on_dx[i] * dx;
	std::cout << dyb[i] << std::endl;
    }

    Polygon<float_> Channel = getNonParallelChannel<float_>(1/(float_)2, dyb[KeyIndices[2]], dyt, H[KeyIndices[1]] - dyt);

    //Polygon<float_> Channel = getNPChannel_ByAngle<float_>(1/(float_)2, dyb[KeyIndices[2]], dyt, H[KeyIndices[1]] - dyt);
    
    Channel.Write("Data/" + Key + "BilliardTable.dat");
    Obstacles[0] = std::make_unique<Polygon<float_>>(Channel);
    Billiard<float_> Table(Obstacles);
    // ======================= BILLIARD DEFINITION END ==========================
   
    // SIMULATION DETAILS
    float_ T = (float_)1e5;
    std::vector<float_> Time = seq<float_>(T, 1000);
    std::vector<float_> TIME = Time;

    Simulation<float_> Sim {
        .Time = TIME,
        .L = {0},
        .TimesToGet = { 1e2, 1e3, 1e4, 1e5 },
        .Symbols = {""},
        .DataDestination = "Data/" + Key,
        .RunTime = ".",
        .MomentStep = (float_)1,
        .TotalTime = TIME.back(),
        .nInitialConditions = 1000,
        .nCollisions = 0,
        .nMoments = 7};


    Writer<float_> w_DX(Sim.DataDestination + "DX.dat", true);
    std::vector<Ensemble<float_>> Moments(Sim.nMoments, Ensemble<float_>(Sim.Time.size()));

    std::vector<size_t> Indices = getIndicesOfVectorValue(Sim.TimesToGet, Sim.Time);
    std::cout << "Collecting displacements at times t = " << index(Sim.Time, Indices) 
	      << " with indices = " << Indices << std::endl;

    for (size_t i = 0; i < (size_t)Sim.nInitialConditions; ++i) {
        if (10 * i % Sim.nInitialConditions == 0 && i != 0)
            std::cout << "Simulated " << i << " initial conditions." << std::endl;

        Particle<float_> P = Table.getInitialState(0, -1, Sim.TotalTime);
        Trajectory<float_> tTraj = Table.getTrajectory(P, Sim.Time);

        std::vector<float_> xDis = tTraj.getDisplacements_x();
	w_DX.WriteRowVector(xDis, Indices);

        std::vector<std::vector<float_>> DXp = tTraj.getDisplacementPowers_x(MomentOrders);
        for (int j = 0; j < Sim.nMoments; ++j)
            Moments[j].Push(DXp[j]);
    }
    
    for (size_t i = 0; i < Moments.size(); ++i) {
	WriteVector<float_>(Sim.DataDestination + "Moment.dat", Moments[i].PullMeans());
    }



    //std::vector<float_> TIME = seq<float_>(2000, 20000, 6);
    //std::vector<float_> ttg = TIME;
    //ttg.erase(ttg.begin());
    //Simulation<float_> Sim {
    //    .Time = TIME,
    //    .L = {0},
    //    .TimesToGet = ttg,
    //    .Symbols = {""},
    //    .DataDestination = "Data/" + Key,
    //    .RunTime = ".",
    //    .MomentStep = (float_)1,
    //    .TotalTime = TIME.back(),
    //    .nInitialConditions = 100000,
    //    .nCollisions = 0,
    //    .nMoments = 7};
    //Sim.Print();
    //getDisplacements_RS(Table, Sim, (float_)0);

    
    //Simulation<float_> Sim {
    //    .Time = TIME,
    //    .L = {0},
    //    .TimesToGet = { 1e2, 1e3, 1e4, 1e5, 1e6, 1e7 },
    //    .Symbols = {""},
    //    .DataDestination = "Data/" + Key,
    //    .RunTime = ".",
    //    .MomentStep = (float_)1,
    //    .TotalTime = TIME.back(),
    //    .nInitialConditions = 1000,
    //    .nCollisions = 0,
    //    .nMoments = 7};
    //Sim.Print();
    //getDisplacementSummary_t(Table, Sim);


    //WriteTrajectory_All<float_>(Table, Sim);
    //VelocityCorrelations<float_>(Table, Sim);  
    //ContinuousSystemSummary<float_>(Table, Sim);  

    //Simulation<float_> FiniteSim {
    //    .Time = Time,
    //    .L = {(float_)100},
    //    .Symbols = {""},
    //    .DataDestination = "Data/" + Key,
    //    .RunTime = ".",
    //    .MomentStep = 1,
    //    .TotalTime = (float_)Time.back(),
    //    .nInitialConditions = 50000,
    //    .nCollisions = 0,
    //    .nMoments = 10};

    //FiniteSim.Write();
    //EscapeToAbsorbingBoundary<float_>(Table, FiniteSim);  
    //FirstPassageTimes<float_>(Table, FiniteSim);
    //ForcedEscapeToAbsorbingBoundary<float_>(Table, FiniteSim);  
    //MeanFirstPassage(Table, FiniteSim);

   return 0;
}

