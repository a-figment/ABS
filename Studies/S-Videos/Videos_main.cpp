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
    DataKeys DK({"t"}, {1}, {"t"});
    std::string Key = DK.MakeKey(KeyIndices);
    DK.Print();
    // ======================= BOILER PLATE END ==========================
    
    // ======================= BILLIARD DEFINITION START ==========================
    
    // Constants
    float_ PHI = boost::math::constants::phi<float_>();
    float_ PI = boost::math::constants::pi<float_>();


    // ==== Billiards ====

    // Parallel channel
    //std::vector<std::unique_ptr<Obstacle<float_>>> Obstacles(1);
    //Polygon<float_> Channel = getParallelChannel2<float_>(0.56, (1 / (float_)2), 0.38);
    //Obstacles[0] = std::make_unique<Polygon<float_>>(Channel);
    //Billiard<float_> Table(Obstacles);

    // Lorentz
    std::vector<std::unique_ptr<Obstacle<float_>>> Obstacles(2);
    Polygon<float_> Square({(float_)0, (float_)1, (float_)1, (float_)0}, {(float_)0, (float_)0, (float_)1, (float_)1}, {1, 3, 0, 2}, true);

    float_ radius = 3 / ((float_)10);
    Circle<float_> Scatterer(radius, Point<float_>((float_)0.5, (float_)0.5));

    Obstacles[0] = std::make_unique<Polygon<float_>>(Square);
    Obstacles[1] = std::make_unique<Circle<float_>>(Scatterer);
    Billiard<float_> Table(Obstacles);
    // ======================= BILLIARD DEFINITION END ==========================
    
    
    //float_ T = (float_)1e4;
    std::vector<float_> TIME = seq<float_>(1/((float_)100), 25,  2500);
    

    Simulation<float_> Sim {
        .Time = TIME,
        .L = {0},
        .TimesToGet = {0},
        .Symbols = {""},
        //.DataDestination = "/fred/oz055/pathtohpc/ABS/S-BallisticFront/" + Key,
        .DataDestination = "Data/" + Key,
        .RunTime = ".",
        .MomentStep = (float_)1,
        .TotalTime = TIME.back(),
        .nInitialConditions = 1000,
        .nCollisions = 0,
        .nMoments = 2};
    //Sim.Write();
    WriteTrajectory_Sample<float_>(Table, Sim);
    //getDisplacements_t<float_>(Table, Sim, {1,2,3,4,5,6});  

    //WriteTrajectory_Sample(Table, Sim);

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

