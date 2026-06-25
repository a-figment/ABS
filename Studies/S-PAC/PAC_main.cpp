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
    DataKeys DK({"p", "a", "d"}, {100, 7, 28}, {"Partitions", "Alpha (RAT) 2pi/3, pi/2, pi/3, ... ", "[0-7] FH, [8-27] IH"});
    std::string Key = DK.MakeKey(KeyIndices);
    DK.Print();
    // ======================= BOILER PLATE END ==========================
    
    // ======================= BILLIARD DEFINITION START ==========================
    std::vector<std::unique_ptr<Obstacle<float_>>> Obstacles(1);
    
    // Internal angle definitions
    float_ PHI = boost::math::constants::phi<float_>();
    float_ PI = boost::math::constants::pi<float_>();

    // 
    std::vector<float_> alpha_rational = { (float_)11*PI/((float_)12), 2 * PI / ((float_)3), PI / ((float_)2), PI /((float_)3), PI / (float_)4, PI/(float_)12};

    std::vector<float_> alpha(alpha_rational.size() * 2);
    for (size_t i = 0; i < alpha.size(); ++i) {
        if (i < (alpha_rational.size()))
            alpha[i] = alpha_rational[i];
        else
            alpha[i] = (alpha_rational[i - alpha_rational.size()] * (PHI - (3 / ((float_)(5)))));
    }
    std::cout << alpha.size() << std::endl;
    
    // 
    float_ dy = (1/((float_)2)) * (cos(alpha[KeyIndices[1]] / (float_)2) / sin(alpha[KeyIndices[1]] / (float_)2)); 

    if (KeyIndices[1] <= alpha_rational.size() - 1)  
        std::cout << "Running RATIONAL alpha = " + std::to_string(alpha[KeyIndices[1]] * 180/PI) << std::endl;
    else
        std::cout << "Running IRRATIONAL alpha = " + std::to_string(alpha[KeyIndices[1]] * 180/PI) << std::endl;



    // 23 ratios
    std::vector<float_> ratio1 = { 30/(float_)100, 40/(float_)100, 50/(float_)100, 60/(float_)100, 70/(float_)100, 80/(float_)100, 90/(float_)100 };
    std::vector<float_> ratio2 = seq<float_>(1, 5, 13);
    //std::vector<float_> ratio2 = seq<float_>(1, 11, 21);
    std::vector<float_> ratio = ratio1;
    ratio.insert(ratio.end(), ratio2.begin(), ratio2.end());
    sort(ratio.begin(), ratio.end());
    //std::vector<float_> porewidths = { dy / 4, dy / 2, 3 * dy/ 4, dy, dy + dy/4, dy + dy/2, dy + 3*dy/4, 2*dy, 3 * dy, 4  * dy, 5 * dy, 6 * dy, 7 * dy, 8 * dy, 9 * dy, 10 * dy };

    //WriteVector<float_>("Data/pw.dat", porewidths);
    //std::vector<float_> rat = { 25/(float_)100, 50/(float_)100, 75/(float_)100, 100/(float_)100, 125/(float_)100, 150/(float_)100, 175/(float_)100, 200/(float_)100, (float_)3, (float_)4, (float_)5, (float_)6, (float_)7, (float_)8, (float_)9, (float_)10 };


    std::vector<float_> porewidths(ratio.size(), 0); //= { - 1/((float_)2), (float_)1 };
    for (int i = 0; i < porewidths.size(); ++i) {
        porewidths[i] = ratio[i] * dy;
    }
    
    //Polygon<float_> Channel = getParallelChannel2<float_>(dy, (1 / (float_)2), porewidths[KeyIndices[2]]);

    // Non-parallel channel
    std::vector<float_> H = { 207 / ((float_)100), 127 / ((float_)100), 117 / ((float_)100), 107 / ((float_)100) };
    float_ dyb = 45/((float_)100);
    float_ dyt = 77/((float_)100);
    Polygon<float_> Channel = getNonParallelChannel<float_>((float_)1, dyb, dyt, H[KeyIndices[1]] - dyt);
    
    Obstacles[0] = std::make_unique<Polygon<float_>>(Channel);
    Billiard<float_> Table(Obstacles);
    // ======================= BILLIARD DEFINITION END ==========================
    
    std::vector<float_> TIME = seq<float_>((float_)200, (float_)500000,  9997);
    
    Simulation<float_> Sim {
        .Time = TIME,
        .L = {0},
        .TimesToGet = { (float_)1000, (float_)10000, (float_)100000, (float_)500000 },
        .Symbols = {""},
        //.DataDestination = "/fred/oz055/jorchard_data/ABS/S-PAC/" + Key,
        .DataDestination = "Data/" + Key,
        .RunTime = ".",
        .MomentStep = 2,
        .TotalTime = TIME.back(),
        .nInitialConditions = 1000,
        .nCollisions = 0,
        .nMoments = 1};
    //Sim.Write();
    
    Sim.Print();
    //getPAC_t<float_>(Table, Sim);  
    getDisplacementSummary_t(Table, Sim);
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

