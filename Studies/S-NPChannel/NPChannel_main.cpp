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
    DataKeys DK({"p", "b"}, {1000, 5}, {"Partitions", "Boundary parameter"});
    std::string Key = DK.MakeKey(KeyIndices);
    DK.Print();
    // ======================= BOILER PLATE END ==========================
    
    // ======================= BILLIARD DEFINITION START ==========================
    std::vector<std::unique_ptr<Obstacle<float_>>> Obstacles(1);
    
    // Internal angle definitions
    //float_ PHI = boost::math::constants::phi<float_>();
    float_ PI = boost::math::constants::pi<float_>();

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

    std::vector<float_> H = { 207 / ((float_)100), 127 / ((float_)100), 117 / ((float_)100), 107 / ((float_)100), 97/(float_)100 };
    float_ dyb = 45/((float_)100);
    float_ dyt = 77/((float_)100);

    //Polygon<float_> Channel = getNonParallelChannel<float_>((float_)1, dyb, dyt, H[KeyIndices[1]] - dyt);
    float_ DYB = 0.2295471036;
    float_ DYT = 0.6295471036;
    float_ d = 0.2;
    float_ alpha = PI - 2*atan(2*DYB);
    float_ beta_t = atan(2*(DYT - d));
    std::cout << beta_t * 180/PI << std::endl;
    float_ normal = PI - beta_t + PI/(float_)2;
    //Polygon<float_> Channel = getNonParallelChannel<float_>(1/(float_)2, DYB, DYT-d, d);

    Polygon<float_> Channel = getNonParallelChannel<float_>(1/(float_)2, dyb, dyt, H[KeyIndices[1]] - dyt);
    
    Obstacles[0] = std::make_unique<Polygon<float_>>(Channel);
    Billiard<float_> Table(Obstacles);


    std::cout << Table.getDiameter() << std::endl;
    // ======================= BILLIARD DEFINITION END ==========================
    
    
    float_ T = 1e3;
    std::vector<float_> Time = seq<float_>(0, T, 1001);
    std::vector<float_> TIME = Time;


    Simulation<float_> Sim {
        .Time = TIME,
        .L = {0},
        .TimesToGet = { (float_)1e3, (float_)1e4, (float_)1e5 },
        .Symbols = {""},
        //.DataDestination = "/fred/oz055/jorchard_data/ABS/S-NPChannel/" + Key,
        .DataDestination = "Data/" + Key,
        .RunTime = ".",
        .MomentStep = 1/(float_)10,
        .TotalTime = TIME.back(),
        .nInitialConditions = 1,
        .nCollisions = 0,
        .nMoments = 35
    };
    

    Particle<float_> P = Table.getInitialState(0, -1, Sim.TotalTime);
    //P.Vel = Point<float_>(1, 0);
    //P.Pos = Point<float_>(-0.5, 0.135);
    float_ x0 = -0.4;
    float_ y0 = 0.3;
    P.Pos = Point<float_>(x0, y0);
    float_ angleToVertex = atan((dyb - y0)/(-x0));
    float_ epsil = PI/(float_)360;
    P.Vel = Point<float_>(cos(angleToVertex+epsil), sin(angleToVertex+epsil));
    Trajectory<float_> tTraj = Table.getTrajectory(P, 30);
    WriteVector(Sim.DataDestination + "x1.dat", tTraj.getLiftedPositions_x());
    WriteVector(Sim.DataDestination + "y1.dat", tTraj.getLiftedPositions_y());

    return 0;
    // Image for thesis
    //Writer<float_> w_X(Sim.DataDestination + "x.dat", true);
    //Writer<float_> w_Y(Sim.DataDestination + "y.dat", true);
    //for (size_t i = 0; i < (size_t)Sim.nInitialConditions; ++i) {
    //    Particle<float_> P = Table.getIC_Rectangle(0, -1, Sim.TotalTime, Point<float_>(-0.4, 0.25), Point<float_>(-0.3, 0.35));
    //    P.Vel = Point<float_>(cos(20*PI/180), sin(20*PI/180));
    //	Trajectory<float_> tTraj = Table.getTrajectory(P, Sim.Time);

    //    w_X.WriteRowVector(tTraj.getLiftedPositions_x());
    //    w_Y.WriteRowVector(tTraj.getLiftedPositions_y());
    //}

    std::vector<float_> MomentOrders = seq<float_>(Sim.MomentStep, Sim.MomentStep * Sim.nMoments, Sim.nMoments);
    std::cout << "Collecting moments at orders p = " << MomentOrders << std::endl;


    Writer<float_> w_DX(Sim.DataDestination + "DX.dat", true);
    std::vector<Ensemble<float_>> Moments(Sim.nMoments, Ensemble<float_>(Sim.Time.size()));

    std::vector<size_t> Indices = getIndicesOfVectorValue(Sim.TimesToGet, Sim.Time);
    std::cout << "Collecting displacements at times t = " << index(Sim.Time, Indices) << " with indices = " << Indices << std::endl;

    for (size_t i = 0; i < (size_t)Sim.nInitialConditions; ++i) {
        if (10 * i % Sim.nInitialConditions == 0 && i != 0)
            std::cout << "Simulated " << i << " initial conditions." << std::endl;

        Particle<float_> P = Table.getInitialState(0, -1, Sim.TotalTime);
        Trajectory<float_> tTraj = Table.getTrajectory(P, Sim.Time);

        std::vector<float_> xDis = tTraj.getDisplacements_x();
	w_DX.WriteRowVector(xDis, Indices);

        //std::vector<std::vector<float_>> DXp = tTraj.getDisplacements_x(Sim.MomentStep, Sim.nMoments);
        std::vector<std::vector<float_>> DXp = tTraj.getDisplacementPowers_x(MomentOrders);
        for (int j = 0; j < Sim.nMoments; ++j)
            Moments[j].Push(DXp[j]);
    }
    

    for (size_t i = 0; i < Moments.size(); ++i)
	WriteVector<float_>(Sim.DataDestination + "Moment.dat", Moments[i].PullMeans());
	//WriteVector<float_>(Sim.DataDestination + "Moment_" + std::to_string(i) + ".dat", Moments[i].PullMeans());

        //Moments[i].Write(Sim.DataDestination + "Moment_" + std::to_string(i) + ".dat");
    //getDisplacementSummary_t(Table, Sim);
    //WriteTrajectory_All<float_>(Table, Sim);
    //VelocityCorrelations<float_>(Table, Sim);  
    //ContinuousSystemSummary<float_>(Table, Sim);  

   return 0;
}

