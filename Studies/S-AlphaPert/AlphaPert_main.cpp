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
typedef float_1000 float_;
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
    DataKeys DK({"p", "a", "d"}, {100, 19, 2}, {"Partitions", "Internal angle: [0-2] Irrational, [3-5] Rational", "Pore widths: [0] Finite horizon, [1] Infinite horizon"});
    std::string Key = DK.MakeKey(KeyIndices);
    DK.Print();
    // ======================= BOILER PLATE END ==========================
    
    // ======================= BILLIARD DEFINITION START ==========================
    std::vector<std::unique_ptr<Obstacle<float_>>> Obstacles(1);
    
    // Internal angle definitions
    float_ PHI = boost::math::constants::phi<float_>();
    float_ PI = boost::math::constants::pi<float_>();
    float_ ETA = (PHI - 3/((float_)5));

    // Contained within [0, pi / 2]
    // std::vector<float_> alpha_irrational = seq((5 * PI) / (2 + 5 * PHI), (5 * PI) / (25 * PHI + 10), 3); // IRRATIONAL
    // std::vector<float_> alpha_rational = seq(PI / (float_)2, PI / (float_)10, 3); // RATIONAL

    // Contained within [0, pi]
    //std::vector<float_> alpha_rational = { 2 * PI / ((float_)3), PI / ((float_)2), PI /((float_)3) };
    //std::vector<float_> alpha_irrational(alpha_rational.size());
    //for (size_t i = 0; i < alpha_irrational.size(); ++i) 
    //    alpha_irrational[i] = (alpha_rational[i] * (PHI - (3 / ((float_)(5)))));

    //// 0-2 IRRATIONAL, 3-6 RATIONAL
    //std::vector<float_> alpha = alpha_irrational;
    //alpha.insert(alpha.end(), alpha_rational.begin(), alpha_rational.end());

    std::vector<float_> porewidths = { - 1/((float_)2), (float_)1 };
    

    //// Epsilon neighbourhood about alpha
    //float_ base = (PI / ((float_)2));
    //std::vector<float_> pwrs = seq<float_>(2, 10, 9);
    //std::vector<float_> alpha;// (2*pwrs.size() + 1);
    //for (size_t i = 0; i < pwrs.size(); ++i)
    //    alpha.push_back(base - pow((float_)10, -pwrs[i]));
    //alpha.push_back(base);
    //for (size_t i = 0; i < pwrs.size(); ++i)
    //    alpha.push_back(base + pow((float_)10, -pwrs[i]));


    ///// ALPHA APPROACHING PI
    //std::vector<float_> alpha = seq<float_>(toRad<float_>(150), toRad<float_>(179), DK.DataKey.second[1]);

    std::vector<float_> alpha = { toRad<float_>(10), toRad<float_>(20), toRad<float_>(45), toRad<float_>(60), toRad<float_>(90), toRad<float_>(115), toRad<float_>(135), toRad<float_>(150) };

    //Polygon<float_> Channel = getParallelChannel<float_>(alpha[KeyIndices[1]] * (PHI - 3/((float_)5)), (1 / (float_)2) , porewidths[KeyIndices[2]]);

    Polygon<float_> Channel = getParallelChannel<float_>(alpha[KeyIndices[1]], (1 / (float_)2) , porewidths[KeyIndices[2]]);

    
    //exit(EXIT_FAILURE);
    Channel.Write("Data/" + Key + "BilliardTable.dat");
    Obstacles[0] = std::make_unique<Polygon<float_>>(Channel);
    Billiard<float_> Table(Obstacles);
    // ======================= BILLIARD DEFINITION END ==========================
    
    //std::vector<float_> Time = getExponentialTimes<float_>(T, 1000),
    
    float_ T = (float_)1e5;
    //std::vector<float_> Time = seq<float_>(1/((float_)10), T,  10000);
    std::vector<float_> Time = getExponentialTimes<float_>(T,  1000);


    Simulation<float_> Sim {
        .Time = Time,
        .L = {0},
        .Symbols = {""},
        .DataDestination = "Data/" + Key,
        .RunTime = ".",
        .MomentStep = (float_)0.5,
        .TotalTime = T,
        .nInitialConditions = 100,
        .nCollisions = 100000,
        .nMoments = 20};
    Sim.Write();
    
    //WriteTrajectory_Collisions<float_>(Table, Sim);
    //VelocityCorrelations<float_>(Table, Sim);  


    Particle<float_> P = Table.getInitialState(0, -1, Sim.Time.back());

    Particle<float_> PC = P;

    P.Vel = Point<float_>(cos(ETA * PI/((float_)4)), sin(ETA * PI/((float_)4)));
    Trajectory<float_> IrratTraj = Table.getTrajectory(PC, Sim.nCollisions); // Collisions only
    IrratTraj.Write(Sim.DataDestination + "IrrationalDirection-" + std::to_string(nDigits) + "-digits.dat");

    Particle<float_> PRat = P;

    P.Vel = Point<float_>(cos(PI/((float_)4)), sin(PI/((float_)4)));
    Trajectory<float_> RatTraj = Table.getTrajectory(PRat, Sim.nCollisions); // Collisions only
    RatTraj.Write(Sim.DataDestination + "RationalDirection-" + std::to_string(nDigits) + "-digits.dat");



   return 0;
}

