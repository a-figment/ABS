#ifndef MISC_INCLUDED
#define MISC_INCLUDED
// Assorted 'billiard simulations' and aglo's with no home

#include <iostream>
#include <vector>
#include <chrono>
#include "Billiard.h"
#include "Simulation.h"
#include "Trajectory.h"
#include "Ensemble.h"
#include "EscapeData.h"
#include "Writer.h"
#include "Helper.h"


template <typename T>
void WriteTrajectory_Time(Billiard<T> &Table, const Simulation<T> &Sim)
{
    uint one_gb_in_bytes = 1000000000;
    if (sizeof(T) * 3 * Sim.Time.size() * Sim.nInitialConditions > one_gb_in_bytes)
        std::cout << "WARNING: Trajectory write estimate is greater than 1GB. " << std::endl;

    int M = Sim.nInitialConditions;
    //if (M > 10)
    //    M = 10;
    for (int i = 0; i < M; ++i)
    {
        Particle<T> P = Table.getInitialState(0, -1, Sim.Time.back());
        Particle<T> PC = P;

        Trajectory<T> ContTraj = Table.getTrajectory(PC, Sim.Time); // Continuous time only
    
	if (Table.Failed)
	    std::cout << "Particle " << i << " failed." << std::endl;

        ContTraj.Write(Sim.DataDestination + "ContTraj_" + std::to_string(i) + ".dat");
    }
}

template <typename T>
void WriteTrajectory_Sample(Billiard<T> &Table, const Simulation<T> &Sim)
{

    uint one_gb_in_bytes = 1000000000;
    if (sizeof(T) * 3 * Sim.Time.size() * Sim.nInitialConditions > one_gb_in_bytes)
        std::cout << "WARNING: Trajectory write estimate is greater than 1GB. " << std::endl;

    T d = Table.Obstacles[0]->getVertices()[5].y;
    T dx = Table.Obstacles[0]->getVertices()[0].x;
    T dy = Table.Obstacles[0]->getVertices()[1].y;

    std::cout << "dy + d/2 = " << dy + d/2 << std::endl;

    // Point<T> Vel(cos(boost::math::constants::pi<T>() / (4 / (T)1)), sin(boost::math::constants::pi<T>() / (4 / (T)1)));
    // for (int i = 0; i < (SX.size() + MX.size() + TX.size()); ++i)
    //
    // For parallel channel
    //std::vector<T> S = seq<T>(0, 1, Sim.nInitialConditions);
    //Beam<T> B(Point<T>(-0.45, d / 4), Point<T>(-0.2, dy + 0.1), 0, 0);
    
    // For Lorentz gas
    std::vector<T> S = seq<T>(0, 1, Sim.nInitialConditions);
    Beam<T> B(Point<T>(0.01, 0.25), Point<T>(0.25, 0.01), 0, 0);

    for (int i = 0; i < Sim.nInitialConditions; ++i)
    {

        Particle<T> P = Table.getInitialState(0, -1, Sim.Time.back());
	P.Pos = B.getPointOnBeam(S[i]);
        //P.Vel = Point<T>(cos(-25 * 3.14/180), sin(-25 * 3.14 / 180));
        P.Vel = Point<T>(cos(43 * 3.14/180), sin(43 * 3.14 / 180));
	
        //Particle<T> PC = P;
	//
        //Particle<T> P = Table.getIC_FixedLine(0, -1, Sim.TotalTime, Point<T>(-0.45, d / 4), Point<T>(-0.2, dy + 0.1));

        //Particle<T> P = Table.getIC_ExpandingBall(0, -1, Sim.TotalTime, Point<T>(dx/2, (-dy/dx)*((dx/2) - dx) + d/2), 0.001);
	
        //Trajectory<T> Traj = Table.getTrajectorySample(P, Sim.Time); // Continuous and distrete time (for plotting)
        //Traj.Write(Sim.DataDestination + "Traj_" + std::to_string(i) + ".dat");

        Trajectory<T> ContTraj = Table.getTrajectory(P, Sim.Time); // Continuous time only
        if (Table.Failed)
            std::cout << "Failed" << std::endl;
        ContTraj.WritePositions(Sim.DataDestination + "Positions_" + std::to_string(i) + ".dat");

        //Trajectory<T> DiscTraj = Table.getTrajectory(PC, Sim.nCollisions); // Collisions only
        //DiscTraj.Write(Sim.DataDestination + "DiscTraj_" + std::to_string(i) + ".dat");
    }
}

template <typename T>
void WriteTrajectory_Collisions(Billiard<T> &Table, const Simulation<T> &Sim)
{
    uint one_gb_in_bytes = 1000000000;
    if (sizeof(T) * 3 * Sim.Time.size() * Sim.nInitialConditions > one_gb_in_bytes)
        std::cout << "WARNING: Trajectory write estimate is greater than 1GB. " << std::endl;

    int M = Sim.nInitialConditions;
    //if (M > 10)
    //    M = 10;
    for (int i = 0; i < M; ++i)
    {
        Particle<T> P = Table.getInitialState(0, -1, Sim.Time.back());
	P.Pos = Point<T>(0, 1.3);
	P.Vel = Point<T>(cos(-40 * 3.14 / 180), sin(-40 * 3.14 / 180));
        Particle<T> PC = P;

        Trajectory<T> DiscTraj = Table.getTrajectory(PC, Sim.nCollisions); // Collisions only
        if (Table.Failed)
            std::cout << "Particle " << i << " failed." << std::endl;
        DiscTraj.Write(Sim.DataDestination + "DiscTraj_" + std::to_string(i) + ".dat");
    }
}


template <typename T>
void WriteTrajectory_Perturb(Billiard<T> &Table, const Simulation<T> &Sim)
{
    uint one_gb_in_bytes = 1000000000;
    if (sizeof(T) * 3 * Sim.Time.size() * Sim.nInitialConditions > one_gb_in_bytes)
        std::cout << "WARNING: Trajectory write estimate is greater than 1GB. " << std::endl;

    //Writer<T> w_DX(Sim.DataDestination + "DX.dat", true);
    
    T d  = Table.Obstacles[0]->getVertices()[5].y;
    T dy = Table.Obstacles[0]->getVertices()[1].y;
    T dx = Table.Obstacles[0]->getVertices()[0].x;
    //T beta = boost::math::constants::pi<T>() - atan2(dy, dx);
    T beta = (T)0;
    T eps = boost::math::constants::pi<T>() / ((T)10000);


    // Sequential distribution in N.H of periodic orbit
    //std::vector<T> ang = seq<T>(beta - eps, beta + eps, Sim.nInitialConditions);
    //ang.push_back(0);
    //std::sort(ang.begin(), ang.end());
    
    // Uniform distribution in N.H of periodic orbit
    std::vector<T> ang(Sim.nInitialConditions, 0);
    Random<T> URAngle(beta - eps, beta + eps);
    for (size_t i = 0; i < Sim.nInitialConditions; ++i) {
        T Theta = URAngle.getUniformRandom();
        ang[i] = Theta;
    }


    T vb = 1 / (1 + (d/dy));
    T thresh = 5/(T)1000;
    std::cout << "dx / t = " << vb << " for the perfect orbit." << std::endl;

    T y0 = d - 1/(T)10;

    if (d > dy) 
	std::cout << "Horizon is infinite." << std::endl;

    std::vector<T> times(Sim.nInitialConditions, 0);
    for (int i = 0; i < Sim.nInitialConditions; ++i)
    {

        if (10 * i % Sim.nInitialConditions == 0 && i != 0)
            std::cout << "Simulated " << i << " initial conditions." << std::endl;

        Particle<T> P = Table.getInitialState(0, -1, Sim.Time.back());
	//P.Pos = Point<T>(dx, y0);
	//P.Vel = Point<T>(cos(ang[i]), sin(ang[i]));
        Particle<T> PC = P;
	Particle<T> PT = P;

	//Trajectory<T> DiscTraj = Table.getTrajectory(PC, Sim.nCollisions); 
        Trajectory<T> ContTraj = Table.getTrajectory(PT, Sim.Time);	   
        if (Table.Failed)
            std::cout << "Particle " << i << " failed." << std::endl;

	std::vector<T> xDis = ContTraj.getDisplacements_x();
	std::vector<T> TSD(xDis.size());
	for (size_t j = 0; j < xDis.size(); ++j) 
	    TSD[j] = xDis[j] / Sim.Time[j];

	//size_t cnt = 10;
	//while (TSD[cnt] > (vb - thresh))
	//    cnt++;
	//times[i] = Sim.Time[cnt];

	WriteVector(Sim.DataDestination + "TSD.dat", TSD);

	//w_DX.WriteRowVector(xDis);
        //ContTraj.Write(Sim.DataDestination + "ContTraj_" + std::to_string(i) + ".dat");
	
	//if (i == 50) {
	//    Trajectory<T> DiscTraj = Table.getTrajectory(PC, Sim.nCollisions); 
	//    DiscTraj.Write(Sim.DataDestination + "Traj-" + std::to_string(i) + ".dat");
	//}
    }
    WriteVector(Sim.DataDestination + "DT.dat", times);
}

// @rand If false, accept the position/velocity arguments
template <typename T>
void WriteOne(Billiard<T> &Table, const Simulation<T> &Sim, const Point<T>& Pos, const Point<T>& Vel, const bool& rand)
{
    uint one_gb_in_bytes = 1000000000;
    if (sizeof(T) * 3 * Sim.Time.size() * Sim.nInitialConditions > one_gb_in_bytes)
        std::cout << "WARNING: Trajectory write estimate is greater than 1GB. " << std::endl;

    
    //T d  = Table.Obstacles[0]->getVertices()[5].y;
    //T dy = Table.Obstacles[0]->getVertices()[1].y;
    //T dx = Table.Obstacles[0]->getVertices()[0].x;
    //T beta = boost::math::constants::pi<T>() - atan2(dy, dx);
    //T eps = boost::math::constants::pi<T>() / (T)720;

    //std::vector<T> ang = seq<T>(beta - eps, beta + eps, Sim.nInitialConditions);


    for (int i = 0; i < Sim.nInitialConditions; ++i)
    {
        Particle<T> P = Table.getInitialState(0, -1, Sim.Time.back());
    
	if (!rand) {
	    P.Pos = Pos;
	    P.Vel = Vel;
	}
        Particle<T> PC = P;
	Particle<T> PT = P;

	//Trajectory<T> DiscTraj = Table.getTrajectory(PC, Sim.nCollisions); 
        Trajectory<T> ContTraj = Table.getTrajectory(PT, Sim.Time);	   


        //DiscTraj.WritePositions(Sim.DataDestination + "DiscTraj-" + std::to_string(i) + ".dat");

	std::vector<T> xDis = ContTraj.getDisplacements_x();
	std::vector<T> TSD(xDis.size());
	for (size_t j = 0; j < xDis.size(); ++j) {
	    TSD[j] = xDis[j] / Sim.Time[j];
	}

        if (Table.Failed)
            std::cout << "Particle " << i << " failed." << std::endl;
        //DiscTraj.Write(Sim.DataDestination + "DiscTraj-" + std::to_string(i) + ".dat");
        //ContTraj.Write(Sim.DataDestination + "ContTraj_" + std::to_string(i) + ".dat");
	WriteVector(Sim.DataDestination + "TSD.dat", TSD);
	

	//WriteVector(Sim.DataDestination + "DXn.dat", DiscTraj.getDisplacements_x());
    }
}



template <typename T>
void WriteTrajectory_All(Billiard<T> &Table, const Simulation<T> &Sim)
{
    uint one_gb_in_bytes = 1000000000;
    if (sizeof(T) * 3 * Sim.Time.size() * Sim.nInitialConditions > one_gb_in_bytes)
        std::cout << "WARNING: Trajectory write estimate is greater than 1GB. " << std::endl;

    int M = Sim.nInitialConditions;
    if (M > 10)
        M = 10;
    for (int i = 0; i < M; ++i)
    {
        Particle<T> P = Table.getInitialState(0, -1, Sim.Time.back());
        Particle<T> PC = P;

        Trajectory<T> Traj = Table.getTrajectorySample(P, Sim.Time); // Continuous and distrete time (for plotting)
        Traj.Write(Sim.DataDestination + "Traj_" + std::to_string(i) + ".dat");

        //Trajectory<T> ContTraj = Table.getTrajectory(PC, Sim.Time); // Continuous time only
        //ContTraj.Write(Sim.DataDestination + "ContTraj_" + std::to_string(i) + ".dat");

        //Trajectory<T> DiscTraj = Table.getTrajectory(PC, Sim.nCollisions); // Collisions only
        //DiscTraj.Write(Sim.DataDestination + "DiscTraj_" + std::to_string(i) + ".dat");
    }
}
// @brief Generates data and basic statistics about a collection of @Sim.nInitialConditions trajectories
//        involving @Sim.nCollisions.
// The following data is written with the @Sim.DataDestination prefix:
//      - MFP.dat - Mean free path for each particle
//      - Moment-i.dat - @Sim.nMoments (displacement) in steps of @Sim.MomentStep
//      - Distribution of trajectory lengths - Considerd at collision numbers 10^i, i = 1,2,...,log10(@Sim.nCollisions)
//      - Velocity autocorrelation -
template <typename T>
inline void CollisionSystemSummary(Billiard<T> &Table, Simulation<T> &Sim)
{
    // Involves all initial conditions
    std::vector<Ensemble<T>> Moments(Sim.nMoments, Ensemble<T>(Sim.nCollisions));
    std::vector<T> MFP(Sim.nInitialConditions, -1);

    // Written per initial condition
    std::vector<T> TrajectoryLengths(static_cast<size_t>(log10(Sim.nCollisions)), -1);
    Ensemble<T> MeanTimeAtCollision(Sim.nCollisions);

    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < Sim.nInitialConditions; ++i)
    {
        if (10 * i % Sim.nInitialConditions == 0 && i != 0)
            std::cout << "Simulated " << i << " initial conditions." << std::endl;

        Particle<T> P0 = Table.getInitialState(0, -1, Sim.TotalTime);
        Trajectory<T> Traj = Table.getTrajectory(P0, Sim.nCollisions);

        // Extract trajectory properties
        MFP[i] = Traj.getMeanFreePath();

        std::vector<T> TrajLengths = Traj.getTimes();
        MeanTimeAtCollision.Push(TrajLengths);

        for (size_t j = 0; j < TrajectoryLengths.size(); ++j)
            TrajectoryLengths[j] = TrajLengths[std::pow(10, (j + 1))];
        WriteVector<T>(Sim.DataDestination + "DS_TrajLengths.dat", TrajectoryLengths);

        std::vector<std::vector<T>> DXp = Traj.getDisplacements_x(Sim.MomentStep, Sim.nMoments);
        for (int j = 0; j < Sim.nMoments; ++j)
            Moments[j].Push(DXp[j]);
    }
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";

    WriteVector<T>(Sim.DataDestination + "DS_MFP.dat", MFP);
    MeanTimeAtCollision.Write(Sim.DataDestination + "DS_MeanTimeAtCollision.dat");
    for (size_t i = 0; i < Moments.size(); ++i)
        Moments[i].Write(Sim.DataDestination + "DS_Moment_" + std::to_string(i) + ".dat");

    Sim.Print(RunTime);
    Sim.Write();
}


// @brief Writes displacements at times t[Indices] 
template <typename T>
inline void getDisplacements_t(Billiard<T> &Table, Simulation<T> &Sim, const std::vector<size_t>& Indices)
{
    std::ofstream Writer_DX(Sim.DataDestination + "DX.dat", std::ios::out | std::ios::app);
    Writer_DX << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;

    //std::ofstream Writer_IC(Sim.DataDestination + "IC.dat", std::ios::out | std::ios::app);

    // Geometry 
    T d  = Table.Obstacles[0]->getVertices()[5].y;
    //T dy = Table.Obstacles[0]->getVertices()[1].y;
    //T dx = Table.Obstacles[0]->getVertices()[0].x;

    Random<T> UR;
    std::vector<T> Heights(Sim.nInitialConditions, -1);
    for (size_t i = 0; i < Sim.nInitialConditions; ++i) 
	Heights[i] = UR.getUniformRandom(1e-12, d - 1e-12);

    //T Vbal = 1 - (1 / (1 + 2*d));
    //std::cout << "Vbal = " << Vbal << std::endl;
    //T eps = 75/(T)1000;

    // Some angles
    //std::vector<T> th = seq<T>(0.2, 2.5 / 4, Sim.nInitialConditions);

    size_t count = 0;
    std::vector<size_t> nCollisions(Sim.nInitialConditions, 0);

    std::cout << "Collecting displacements at times t = ";
    for (size_t i = 0; i < Indices.size(); ++i)
        std::cout << Sim.Time[Indices[i]] << ", ";
    std::cout << std::endl;

    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < Sim.nInitialConditions; ++i)
    {
        if (10 * i % Sim.nInitialConditions == 0 && i != 0)
            std::cout << "Simulated " << i << " initial conditions." << std::endl;

        // ON P.B.
        //Particle<T> P0 = Table.getIC_FixedLine(0, -1, Sim.TotalTime, Point<T>(-dx, (T)0), Point<T>(-dx, d));
        //P0.Vel = Point<T>(1,0);
        //
        //Particle<T> P0 = Table.getIC_FixedPosition(0, -1, Sim.TotalTime, Point<T>(0, dy + 0.0001));
        //P0.Vel = Point<T>(cos(th[i]), sin(th[i]));
        //std::cout << P0 << std::endl;
        
        // Standard distribution
        Particle<T> P0 = Table.getInitialState(0, -1, Sim.TotalTime);
	//P0.Pos.x = -(1/(T)2) + 0.000000001;
	//P0.Pos.y = Heights[i];
        Particle<T> Pn = P0;

        Trajectory<T> Traj = Table.getTrajectory(Pn, Sim.Time);
        //Trajectory<T> Traj = Table.getTrajectory(P0, Sim.nCollisions);
        nCollisions[i] = Traj.nCollisions;

        //
        std::vector<int> flag(Indices.size(), 0);

        // Set and write displacements
        std::vector<T> xDis = Traj.getDisplacements_x();
	bool InThresh = false;
        for (size_t j = 0; j < Indices.size(); ++j) {
	    Writer_DX << xDis[Indices[j]] << " ";
            //T v = abs(xDis[Indices[j]]) / Sim.Time[Indices[j]]; // Make sure v = .5 at all times
            //if (v < 0.57 || v > 0.59)  
    	    InThresh = false;
        }
	Writer_DX << std::endl;
        
        // Saturated at all values
        //if (InThresh && count < 10)
        //{
	//    WriteVector<T>({P0.Pos.x, P0.Pos.y, atan2(P0.Vel.y, P0.Vel.x)}, Writer_IC);
	//    Particle<T> Pc = P0, Pt = P0;
        //    Trajectory<T> Traj1 = Table.getTrajectory(Pc, Sim.nCollisions);
	//    Trajectory<T> Traj2 = Table.getTrajectory(Pt, Sim.Time);
        //    Traj1.WritePositions(Sim.DataDestination + "ColTraj-" + std::to_string(i) + ".dat");
        //    Traj2.WritePositions(Sim.DataDestination + "ConTraj-" + std::to_string(i) + ".dat");
	//    count++;
        //}
    }

    //WriteVector(Sim.DataDestination + "nCollisions.dat", nCollisions);

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";
    std::cout << RunTime << std::endl;
    //Sim.Print(RunTime);
}

// @brief Writes moments of displacement
//template <typename T>
//inline void getMoments_t(Billiard<T> &Table, Simulation<T> &Sim, const T& eps)
//{
//    // Involves all initial conditions
//    std::vector<Ensemble<T>> Moments(Sim.nMoments, Ensemble<T>(Sim.Time.size()));
//
//    std::vector<T> thetas = generateICs(Table, Sim.Time, Sim.nCollisions, Sim.nInitialConditions, Point<T>(100, -1));
//
//    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
//    for (int i = 0; i < Sim.nInitialConditions; ++i)
//    //int i = 0;
//    //while (i < Sim.nInitialConditions)
//    {
//        if (10 * i % Sim.nInitialConditions == 0 && i != 0)
//            std::cout << "Simulated " << i << " initial conditions." << std::endl;
//
//        Particle<T> P0 = Table.getInitialState(0, -1, Sim.TotalTime);
//	P0.Vel = Point<T>(cos(thetas[i], sin(thetas[i])));
//	
//	//T th = std::abs(atan(P0.Vel.y / P0.Vel.x));
//	//T p2 = boost::math::constants::pi<T>()/(T)2;
//	//T p2 = (T)1;
//	//T p3 = 58/(T)100;
//
//
//	//if ((th < p2 && th > p2 - eps) || (th > -p2 && th < -p2+eps) || (th < eps && th > -eps)) {
//	//if ((th < p2 + eps && th > p2 - eps) || (th < p3 - eps && th > p3 + eps)) {
//	//    continue;
//	//} else {
//	    Trajectory<T> Traj = Table.getTrajectory(P0, Sim.Time);
//
//            // Extract trajectory properties
//            std::vector<std::vector<T>> DXp = Traj.getDisplacements_x(Sim.MomentStep, Sim.nMoments);
//            for (int j = 0; j < Sim.nMoments; ++j)
//                Moments[j].Push(DXp[j]);
//	    //i++;
//	//}
//	    
//    }
//    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
//    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";
//
//    for (size_t i = 0; i < Moments.size(); ++i)
//        Moments[i].Write(Sim.DataDestination + "Moment_" + std::to_string(i) + ".dat");
//    //Sim.Print(RunTime);
//}

// @brief Gets continuous time position-position autocorrelations
template <typename T>
inline void getPAC_t(Billiard<T> &Table, Simulation<T> &Sim)
{
    std::vector<size_t> h = { 0, 10, 50, 100, 500, 1000, 5000 }; // Assumes time step of unity

    // Involves all initial conditions
    std::vector<Ensemble<T>> PAC(h.size(), Ensemble<T>(Sim.Time.size() - h[0]));

    std::cout << "Collecting PAC with steps h = ";
    for (size_t i = 0; i < h.size()-1; ++i)
        std::cout << h[i] << ", ";
    std::cout << h.back() << "." << std::endl;

    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < Sim.nInitialConditions; ++i)
    {
        if (10 * i % Sim.nInitialConditions == 0 && i != 0)
            std::cout << "Simulated " << i << " initial conditions." << std::endl;

        Particle<T> P0 = Table.getInitialState(0, -1, Sim.TotalTime);
        Trajectory<T> Traj = Table.getTrajectory(P0, Sim.Time);

        // Get 
        std::vector<T> xDis = Traj.getDisplacements_x();
        for (size_t i = 0; i < h.size(); ++i) {
            std::vector<T> pac_temp(xDis.size() - h[i]); // trailing zeros
            for (size_t j = 0; j < xDis.size() - h[i]; ++j) {
                pac_temp[j] = xDis[j] * xDis[j + h[i]];
            }
            PAC[i].Push(pac_temp);
        }
    }

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";

    for (size_t i = 0; i < PAC.size(); ++i)
        PAC[i].Write(Sim.DataDestination + "PAC_" + std::to_string(i) + ".dat");
    Sim.Print(RunTime);
}

// @brief Looks at all quantifiers of the displacements of particles in continuous time units:
// 1. Displacements at selected times (used for constructing density estimates)
// 2. Moments of displacment 
// 3. Position auto-correlation
// 4. TAMSD
// 5. EA TAMSD
template <typename T>
inline void getDisplacementSummary_t(Billiard<T> &Table, Simulation<T> &Sim)
{
    // Streams 
    Writer<T> w_DX(Sim.DataDestination + "DX.dat", true);
    //Writer<T> wTA(Sim.DataDestination + "TA.dat", true);
    Writer<T> w_TAMSD(Sim.DataDestination + "TAMSD.dat", true);

    Writer<T> w_tTAMSD0(Sim.DataDestination + "tTAMSD-0.dat", true);
    Writer<T> w_tTAMSD1(Sim.DataDestination + "tTAMSD-1.dat", true);
    Writer<T> w_tTAMSD2(Sim.DataDestination + "tTAMSD-2.dat", true);
    Writer<T> w_tTAMSD3(Sim.DataDestination + "tTAMSD-3.dat", true);


    // Moments of displacement
    std::vector<Ensemble<T>> Moments(Sim.nMoments, Ensemble<T>(Sim.Time.size()));

    // 
    std::vector<size_t> szs = {0, 9, 99, 999};		// Specific Delta indices to consider (indexes h)
    std::vector<size_t> h = seq<size_t>(1, 1000, 1000); // * dt = Delta. Must satisfy h.back() < Sim.Time.back()

    // Positiona auto-correlation for lag times h[szs]
    std::vector<Ensemble<T>> PAC(szs.size(), Ensemble<T>(Sim.Time.size() - h.back()));

    // EA TAMSD as a function of \Delta (h)
    Ensemble<T> EATAMSD(h.size()); // Ensemble averaged time averaged mean squared displacement
    
    // 
    std::vector<Ensemble<T>> tEATAMSD(szs.size(),  Ensemble<T>(Sim.Time.size() - h[szs.back()]));     
    std::vector<Ensemble<T>> tEATAMSD2(szs.size(), Ensemble<T>(Sim.Time.size() - h[szs.back()])); 

    // TAMSD for 
    //std::vector<std::vector<T>> tTAMSD(szs.size(), std::vector<T>(Indices.size()));

    std::vector<size_t> Indices = getIndicesOfVectorValue(Sim.TimesToGet, Sim.Time);
    std::cout << "Collecting displacements at times t = ";
    for (size_t i = 0; i < Indices.size(); ++i)
        std::cout << Sim.Time[Indices[i]] << " ";
    std::cout << std::endl;

    T dt = Sim.Time[2] - Sim.Time[1];
    if ((Sim.Time[42] - Sim.Time[41]) != dt)
	std::cout << "WARNING: Observation times are not linear. Correlations and TA's will be affected." << std::endl; 
    else
        std::cout << "Linear time in steps of dt = " << dt << "." << std::endl; 

    //if ((Sim.Time.size() - h[szs.back()]) > Indices.back())
    if (Indices.size() < 3)
	std::cout << "WARNING: At least 3 observation times are needed to write tTAMSD data." << std::endl; 

    //int nFound = 0;
    //T eps = 0.075;
    //T Vbal = 0.4;
    //std::cout << "Finding trajectories with " << Vbal - eps << " < Vbal < " << Vbal + eps << std::endl;


    // Main loop
    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    for (int n = 0; n < Sim.nInitialConditions; ++n)
    {
        if (10 * n % Sim.nInitialConditions == 0 && n != 0)
            std::cout << "Simulated " << n << " initial conditions." << std::endl;

        Particle<T> P0 = Table.getInitialState(0, -1, Sim.TotalTime);
        Trajectory<T> Traj = Table.getTrajectory(P0, Sim.Time);

        // GET MOMENTS OF DISPLACEMENT
	std::vector<std::vector<T>> DXp = Traj.getDisplacements_x(Sim.MomentStep, Sim.nMoments);
        for (int j = 0; j < Sim.nMoments; ++j)
            Moments[j].Push(DXp[j]);
	
	// WRITE DISPLACEMENTS
        std::vector<T> xDis = Traj.getDisplacements_x();
	w_DX.WriteRowVector(xDis, Indices);

	//bool ThreshCond = true;
        //for (size_t j = 0; j < Indices.size(); ++j) {
        //    T v = abs(xDis[Indices[j]]) / Sim.Time[Indices[j]]; // Make sure v = .5 at all times
	//    if (v < Vbal - eps || v > Vbal + eps) { // Only one instance of dx/t < .2 removes the sample
	//	//std::cout << v << std::endl;
	//	ThreshCond = false;
	//	break;
	//    }
	//}
	//if (ThreshCond) {
	//    nFound++;
	//    std::cout << "FOUND " << nFound << std::endl;

	//    std::vector<T> delta2_bar(h.size(),0); // For the current trajectory, average all 
        //    for (size_t i = 0; i < h.size(); ++i) {
	//        T delta2 = 0;
        //        for (size_t j = 0; j < xDis.size() - h[i]; ++j) {
	//	    delta2 += (xDis[j + h[i]] - xDis[j]) * (xDis[j + h[i]] - xDis[j]);
        //        }
	//        delta2_bar[i] = delta2 / (Sim.Time.back() - h[i] * (dt));
        //    }
	//    wTA.WriteVector(delta2_bar);
	////    for (size_t z = 0; z < delta2_bar.size(); ++z)	
	////	Writer_TA << delta2_bar[z] << " ";
	////    Writer_TA << std::endl;

	//    // Set and write displacements
        //    for (size_t j = 0; j < Indices.size(); ++j) 
        //        Writer_DX << xDis[Indices[j]] << " ";
        //    Writer_DX << std::endl;
	//} else {
	//    Sim.nInitialConditions++;
	//}

	for (size_t i = 0; i < szs.size(); ++i) {
            std::vector<T> pac(xDis.size() - h[szs[i]], 0); // trailing zeros
	    for (size_t j = 0; j < xDis.size() - h[szs[i]]; j++) 
                pac[j] = xDis[j] * xDis[j + h[szs[i]]];
            PAC[i].Push(pac);
	}

        // ================ \bar{\delta^2}(\Delta; T) - TAMSD FOR A FIXED TOTAL TIME ==================
	std::vector<T> delta2_bar(h.size(),0);	       // For the current trajectory, average all 
        for (size_t i = 0; i < h.size(); ++i) {
	    T delta2 = 0;
            for (size_t j = 0; j < xDis.size() - h[i]; ++j) {
		delta2 += ((xDis[j + h[i]] - xDis[j]) * (xDis[j + h[i]] - xDis[j]));
            }
	    T time = Sim.Time[Sim.Time.size() - h[i] - 1]; // = Sim.Time.back() - h[i] * dt;
	    //T time = Sim.Time.back() - h[i] * dt;
	    delta2_bar[i] = delta2 / time;
        }
	w_TAMSD.WriteRowVector(delta2_bar, szs);
	EATAMSD.Push(delta2_bar);

        // ============================================================================================


	//// =========== \bar{\delta^2}(t; \Delta) - TAMSD FOR A FIXED \Delta = h*dt =============
	std::vector<std::vector<T>> tTAMSD(szs.size(), std::vector<T>(3, 0));
	for (size_t i = 0; i < szs.size(); ++i) { // Loop over a few values of \Delta = h[szs[i]]

	    std::vector<T> delta2_bar(tEATAMSD[i].nElems, 0); // Truncated by largest delta value
	    std::vector<T> delta2_bar2(tEATAMSD2[i].nElems, 0); // Truncated by largest delta value
	    T delta2 = 0;
	    for (size_t j = h[szs[i]]-1; j < delta2_bar.size(); ++j) {

		T time = Sim.Time[j - (h[szs[i]] - 1)]; // = Sim.Time.back() - h[i] * dt;
		delta2 += ((xDis[j + h[szs[i]]] - xDis[j])*(xDis[j + h[szs[i]]] - xDis[j]));
		
		delta2_bar[j]  = delta2 / time;
		delta2_bar2[j] = delta2_bar[j] * delta2_bar[j];
	    }

	    // Cumulative time average for a fixed Delta
	    //size_t y = xDis.size() - h[szs[i]];
	    //std::vector<T> CumSum(y,0);
	    //std::vector<T> CumSumSq(y,0);

	    //// 
	    //size_t start_idx = 0;
	    //while (Sim.Time[start_idx] < h[szs[i]] * dt) 
	    //    start_idx++;

	    //for (size_t j = start_idx+1; j < y - 1; ++j) { // For t - \Delta > 0, find cumulative time average

	    //    T delta2_bar = (xDis[j + h[szs[i]]] - xDis[j]) * (xDis[j + h[szs[i]]] - xDis[j]);
	    //    T curr_obs_time = Sim.Time[j] - h[szs[i]] * dt;

	    //    // Pre-padded with zeros
	    //    CumSum[j] = CumSum[j-1] + ( delta2_bar / curr_obs_time );
	    //    CumSumSq[j] = CumSum[j] * CumSum[j];
	    //}

	    tTAMSD[i] = { delta2_bar[Indices[0]], delta2_bar[Indices[1]], delta2_bar[Indices[2]] };

	    tEATAMSD[i].Push(delta2_bar);
	    tEATAMSD2[i].Push(delta2_bar2);
	}
	w_tTAMSD0.WriteRowVector(tTAMSD[0]);
	w_tTAMSD1.WriteRowVector(tTAMSD[1]);
	w_tTAMSD2.WriteRowVector(tTAMSD[2]);
	w_tTAMSD3.WriteRowVector(tTAMSD[3]);
	// ======================================================================

    }

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";

    for (size_t i = 0; i < Moments.size(); ++i)
        Moments[i].Write(Sim.DataDestination + "Moment_" + std::to_string(i) + ".dat");

    for (size_t i = 0; i < tEATAMSD.size(); ++i) {
	std::vector<T> _means = tEATAMSD[i].PullMeans();
	std::vector<T> means2 = tEATAMSD2[i].PullMeans();
	std::vector<T> EB(_means.size(), 0);
	for (size_t j = 0; j < EB.size(); ++j) {
	    EB[j] = (means2[j] - _means[j]*_means[j]) / (_means[j] * _means[j]);
	}
	WriteVector<T>(Sim.DataDestination + "EB.dat", EB);
        tEATAMSD[i].Write(Sim.DataDestination + "tEATAMSD-" + std::to_string(i) + ".dat");
        tEATAMSD2[i].Write(Sim.DataDestination + "tEATAMSD2-" + std::to_string(i) + ".dat");
    }

    for (size_t i = 0; i < szs.size(); ++i) 
        PAC[i].Write(Sim.DataDestination + "PAC_" + std::to_string(i) + ".dat");
    EATAMSD.Write(Sim.DataDestination + "EATAMSD.dat");

    Sim.Print(RunTime);
}


template <typename T>
inline void getVelocities_t(Billiard<T> &Table, Simulation<T> &Sim)
{
    Writer<T> _Theta(Sim.DataDestination + "Theta.dat", true);
    //T d = Table.Obstacles[0]->getVertices()[5].y;
    //T dx = Table.Obstacles[0]->getVertices()[0].x;

    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < Sim.nInitialConditions; ++i)
    {
        if (10 * i % Sim.nInitialConditions == 0 && i != 0)
            std::cout << "Simulated " << i << " initial conditions." << std::endl;

        //Particle<T> P0 = Table.getIC_FixedLine(0, -1, Sim.TotalTime, Point<T>(dx, (T)0), Point<T>(dx, d));
        Particle<T> P0 = Table.getInitialState(0, -1, Sim.TotalTime);
        //
        Trajectory<T> Traj = Table.getTrajectory(P0, Sim.Time);

        // Set and write displacements
	_Theta.WriteRowVector(Traj.getVelocities_a());
    }

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";
    Sim.Print(RunTime);
}

// @brief Gets the flight times of billiard segments which pass 1 or more E.C's without a collision
template <typename T>
inline void getIHFlightTimes(Billiard<T> &Table, Simulation<T> &Sim)
{
    //Writer<T> _IHFT(Sim.DataDestination + "IHFT.dat", true);
    Writer<T> _X(Sim.DataDestination + "X.dat", true);
    //T d = Table.Obstacles[0]->getVertices()[5].y;
    //T dx = Table.Obstacles[0]->getVertices()[0].x;

    //std::vector<T> FT;
    //FT.reserve(10000000);
    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();

    T lb = 20;
    uint i = 0;
    while (i < Sim.nInitialConditions) {
        if (10 * i % Sim.nInitialConditions == 0 && i != 0)
            std::cout << "Simulated " << i << " initial conditions." << std::endl;

        Particle<T> P0 = Table.getInitialState(0, -1, Sim.TotalTime);
	Particle<T> P1 = P0;
        //Trajectory<T> Traj = Table.getTrajectory(P0, Sim.nCollisions);

	Trajectory<T> cTraj = Table.getTrajectory(P1, Sim.Time);
	_X.WriteRowVector(cTraj.getLiftedPositions_x());
	//_X.WriteRowVector(Traj.getLiftedPositions_x());

        // Set and write displacements
	//
	//std::vector<T> lfp = Traj.getLiftedFreePaths();
	//_IHFT.WriteRowVector(lfp);

//	bool can_consider = false;
//	for (size_t j = 0; j < lfp.size(); ++j) {
//	    if (lfp[j] > lb) {
//		can_consider = true;
//		break;
//		//FT.push_back(lfp[j]);
//		//i++;
//	    }
//	}
//
//	if (can_consider) {
//	    //Trajectory<T> cTraj = Table.getTrajectory(P1, Sim.Time);
//	    //_X.WriteRowVector(cTraj.getLiftedPositions_x());
//
//	    _X.WriteRowVector(Traj.getLiftedPositions_x());
//	    _IHFT.WriteRowVector(lfp);
//	    i++;
//	}
        i++;

    }

    //_IHFT.WriteColVector(FT);

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";
    std::cout << RunTime << std::endl;
}


// @brief Gets the slope of the segment joining two collisions. 
template <typename T>
inline void getCollisionVelocity(Billiard<T> &Table, Simulation<T> &Sim)
{
    Writer<T> _M(Sim.DataDestination + "M.dat", true);
    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < Sim.nInitialConditions; ++i) {
        if (10 * i % Sim.nInitialConditions == 0 && i != 0)
            std::cout << "Simulated " << i << " initial conditions." << std::endl;

        Particle<T> P0 = Table.getInitialState(0, -1, Sim.TotalTime);
        Trajectory<T> Traj = Table.getTrajectory(P0, Sim.nCollisions);

	std::vector<Point<T>> Pos = Traj.getLiftedPositions();
	std::vector<T> M(Pos.size()-1, 0);
	for (size_t j = 0; j < Pos.size()-1; ++j) {
	    M[j] = ((Pos[j+1].y - Pos[j].y) / (Pos[j+1].x - Pos[j].x));
	}

	_M.WriteRowVector(M);
    }


    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";
    std::cout << RunTime << std::endl;
}




// @brief All velocities are initially horizontal with positions that span the vertical coordinates of the bounding box
template <typename T>
inline void getHPO(Billiard<T> &Table, Simulation<T> &Sim)
{
    Writer<T> _LFP(Sim.DataDestination + "LFP.dat", true);
    Writer<T> _DX(Sim.DataDestination + "DX.dat", true);
    Writer<T> _POS(Sim.DataDestination + "IPOS.dat", true);
    Writer<T> _TH(Sim.DataDestination + "TH.dat", true);

    // Geometry 
    T d  = Table.Obstacles[0]->getVertices()[5].y;
    T dy = Table.Obstacles[0]->getVertices()[1].y;
    T dx = Table.Obstacles[0]->getVertices()[0].x;
    T m  = dy / (-dx);

    // Horizontal velocity. By symmetry we don't need the other direction
    Point<T> Vel(1, 0);
    int n = std::floor(Sim.nInitialConditions / 2); 
    int factor = Sim.nInitialConditions;

    // Define two beams corresponding to channel segments
    T ub = dy;
    if (d < dy)
	ub = d;
    Beam<T> B1(Point<T>(dx,0), Point<T>(dx, ub), 0, 0);
    int nPointOnB1 = std::floor(factor * B1.Length);

    Beam<T> B2(Point<T>(dx,d), Point<T>(0, dy + d), 0, 0);
    int nPointOnB2 = std::floor(factor * B2.Length);

    std::vector<Particle<T>> Particles(nPointOnB1 + nPointOnB2 - 1);
    for (size_t i = 1; i < Particles.size(); ++i) {
	if (i < nPointOnB1)
	    Particles[i] = Particle<T>(Point<T>(B1.getPointOnBeam( (T)i / nPointOnB1 )), Vel, 0, -1, Sim.Time.back());
	else
	    Particles[i] = Particle<T>(Point<T>(B2.getPointOnBeam( (T)(i - nPointOnB1 + 1) / nPointOnB2 )), Vel, 0, -1, Sim.Time.back());
	_POS.WriteRowVector({ Particles[i].Pos.x, Particles[i].Pos.y });
    }
    
    std::vector<size_t> Indices = getIndicesOfVectorValue(Sim.TimesToGet, Sim.Time);
    std::cout << "Collecting displacements at times t = " << index(Sim.Time, Indices) << " with indices = " << Indices << std::endl;

    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    for (int i = 1; i < Particles.size(); ++i)
    {
        if (10 * i % Particles.size() == 0 && i != 0)
            std::cout << "Simulated " << i << " initial conditions." << std::endl;

        Particle<T> P0 = Particles[i];
        Trajectory<T> Traj = Table.getTrajectory(P0, Sim.nCollisions);
        //Trajectory<T> Traj = Table.getTrajectory(P0, Sim.Time);
	

	// Angle + length of segments over one period should identify the orbit uniquely		
	// Although.. length alone is probably enough
	std::vector<T> LFP = Traj.getLiftedFreePaths();
	_LFP.WriteRowVector(LFP);

        Traj.WritePositions(Sim.DataDestination + "ColTraj-" + std::to_string(i) + ".dat");

	std::vector<T> a = Traj.getVelocities_a();
	_TH.WriteRowVector(a);
        //std::vector<T> xDis = Traj.getDisplacements_x();
	//_DX.WriteRowVector(xDis, Indices);
        //for (size_t i = 0; i < h.size(); ++i) {
        //    std::vector<T> pac_temp(xDis.size() - h[i]); // trailing zeros
        //    for (size_t j = 0; j < xDis.size() - h[i]; ++j) {
        //        pac_temp[j] = xDis[j] * xDis[j + h[i]];
        //    }
        //    PAC[i].Push(pac_temp);
        //}
    }

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";

    std::cout << RunTime << std::endl;
}


// @brief All velocities are perturbed from an initially horizontal directoin with positions that are random
template <typename T>
inline void getPPO(Billiard<T> &Table, Simulation<T> &Sim)
{
    std::ofstream Writer_DX(Sim.DataDestination + "DX.dat", std::ios::out | std::ios::app);
    Writer_DX << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;

    std::ofstream Writer_POS(Sim.DataDestination + "POS.dat", std::ios::out | std::ios::app);
    Writer_POS << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;

    // Geometry 
    T d  = Table.Obstacles[0]->getVertices()[5].y;
    T dy = Table.Obstacles[0]->getVertices()[1].y;
    T dx = Table.Obstacles[0]->getVertices()[0].x;
    T m  = dy / (-dx);

    // Horizontal velocity. By symmetry we don't need the other direction
    Point<T> Vel(1, 0);
    
    int n = std::floor(Sim.nInitialConditions / 2); 
    //std::cout << "n = " << n << std::endl;

    int factor = Sim.nInitialConditions;

    // Define two beams corresponding to channel segments

    Beam<T> B5(Point<T>(dx,0), Point<T>(dx, d), 0, 0);
    int nPointOnB5 = std::floor(factor * B5.Length);

    Beam<T> B4(Point<T>(dx,d), Point<T>(0, dy + d), 0, 0);
    int nPointOnB4 = std::floor(factor * B4.Length);

    std::vector<Particle<T>> Particles(nPointOnB4 + nPointOnB5 - 1);
    for (size_t i = 1; i < Particles.size(); ++i) {
	if (i < nPointOnB5)
	    Particles[i] = Particle<T>(Point<T>(B5.getPointOnBeam( (T)i / nPointOnB5 )), Vel, 0, -1, Sim.Time.back());
	else
	    Particles[i] = Particle<T>(Point<T>(B4.getPointOnBeam( (T)(i - nPointOnB5 + 1) / nPointOnB4 )), Vel, 0, -1, Sim.Time.back());
	Writer_POS << Particles[i].Pos.x << " " << Particles[i].Pos.y << std::endl;	
    }
    
    std::vector<size_t> Indices = getIndicesOfVectorValue(Sim.TimesToGet, Sim.Time);

    std::cout << "Indices found: ";
    for (size_t i = 0; i < Indices.size(); ++i)
        std::cout << Indices[i] << " ";
    std::cout << std::endl;
    std::cout << "Collecting displacements at times t = ";
    for (size_t i = 0; i < Indices.size(); ++i)
        std::cout << Sim.Time[Indices[i]] << " ";
    std::cout << std::endl;

    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    for (int i = 1; i < Particles.size(); ++i)
    {
        if (10 * i % Sim.nInitialConditions == 0 && i != 0)
            std::cout << "Simulated " << i << " initial conditions." << std::endl;

        Particle<T> P0 = Particles[i];
        Trajectory<T> Traj = Table.getTrajectory(P0, Sim.nCollisions);
        //Trajectory<T> Traj = Table.getTrajectory(P0, Sim.Time);
	
        Traj.WritePositions(Sim.DataDestination + "ColTraj-" + std::to_string(i) + ".dat");

        // Get 
        //std::vector<T> xDis = Traj.getDisplacements_x();
        //for (size_t j = 0; j < xDis.size(); ++j)
        //    Writer_DX << xDis[Indices[j]] << " ";
        //Writer_DX << std::endl;
        //for (size_t i = 0; i < h.size(); ++i) {
        //    std::vector<T> pac_temp(xDis.size() - h[i]); // trailing zeros
        //    for (size_t j = 0; j < xDis.size() - h[i]; ++j) {
        //        pac_temp[j] = xDis[j] * xDis[j + h[i]];
        //    }
        //    PAC[i].Push(pac_temp);
        //}
    }

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";

    Sim.Print(RunTime);
}



// e.g. for Lorentz / Ehrenfest
template <typename T>
inline void ContinuousSystemSummary2D(Billiard<T> &Table, Simulation<T> &Sim)
{
    std::ofstream Writer_R(Sim.DataDestination + "R.dat", std::ios::out | std::ios::app);
    Writer_R << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;

    // Involves all initial conditions
    std::vector<Ensemble<T>> Moments(Sim.nMoments, Ensemble<T>(Sim.Time.size()));

    //Ensemble<T> MeanPartialSum(Sim.Time.size());
    //std::vector<int> Indices = {1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000};
    //
    std::vector<int> Indices = {401, 602, 803, 1003};

    std::cout << "Collecting displacements at times t = " << index(Sim.Time, Indices) << " with indices = " << Indices << std::endl;

    std::vector<T> xDisplacement(Indices.size(), -1);

    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < Sim.nInitialConditions; ++i)
    {
        if (10 * i % Sim.nInitialConditions == 0 && i != 0)
            std::cout << "Simulated " << i << " initial conditions." << std::endl;

        Particle<T> P0 = Table.getInitialState(0, -1, Sim.TotalTime);
        Trajectory<T> Traj = Table.getTrajectory(P0, Sim.Time);

        // Extract trajectory properties
        std::vector<std::vector<T>> Rp = Traj.getDisplacements_xy(Sim.MomentStep, Sim.nMoments);
        for (int j = 0; j < Sim.nMoments; ++j)
            Moments[j].Push(Rp[j]);

        std::vector<T> xDis = Traj.getDisplacements_xy();

        // Set and write displacements
        for (size_t j = 0; j < xDisplacement.size(); ++j)
            xDisplacement[j] = xDis[Indices[j] - 1];
        for (size_t i = 0; i < xDisplacement.size(); ++i)
            Writer_R << xDisplacement[i] << " ";
        Writer_R << std::endl;

        //WriteVector<T>(Sim.DataDestination + "xPositions.dat", xPositions);
        //WriteVector<T>(Sim.DataDestination + "xDisplacements.dat", xDisplacement);
        //Traj.MakeParticle(0).WriteState(Sim.DataDestination + "InitialConditions.dat");
    }

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";

    for (size_t i = 0; i < Moments.size(); ++i)
        Moments[i].Write(Sim.DataDestination + "Moment_" + std::to_string(i) + ".dat");
    Sim.Print(RunTime);
    Sim.Write();
}


template <typename T>
void EvolveBeam(Beam<T>& B, Billiard<T>& Table, Simulation<T> &Sim, const size_t& nPointsOnBeam, const size_t& nthBeam)
{
    Writer<int> w_ID(Sim.DataDestination + "ID.dat", true);
    Writer<long int> w_CID(Sim.DataDestination + "CID.dat", true);

    Point<T> Vel(cos(B.VelAngle), sin(B.VelAngle));
    for (size_t i = 1; i < nPointsOnBeam; ++i)
    {
        T t = (T)i / (T)nPointsOnBeam;

        Point<T> Pos = B.getPointOnBeam(t);
        Particle<T> P(Pos, Vel, 0, -1, Sim.TotalTime);

        Trajectory<T> Traj = Table.getTrajectory(P, Sim.nCollisions);

        std::vector<int> IDs = Traj.getSegmentIDs();
        std::vector<long int> CIDs = Traj.getCellIDs_x();
	w_ID.WriteRowVector(IDs);
	w_CID.WriteRowVector(CIDs);

        //std::vector<Point<T>> Poss = Traj.getLiftedPositions();
        //for (size_t i = 0; i < Poss.size(); ++i)
        //    Writer_R << Poss[i].x << " " << Poss[i].y << std::endl;
    }
}


// @brief Rejection samples displacements based on ballistic speed
template <typename T>
inline void getDisplacements_RS(Billiard<T> &Table, Simulation<T> &Sim, const T& Vbal)
{
    std::ofstream Writer_DX(Sim.DataDestination + "DX.dat", std::ios::out | std::ios::app);
    Writer_DX << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;

    // Geometry 
    //T d  = Table.Obstacles[0]->getVertices()[5].y;
    //T dy = Table.Obstacles[0]->getVertices()[1].y;
    //T dx = Table.Obstacles[0]->getVertices()[0].x;

    //T Vbal = 2/(T)10;

    // Some angles
    //std::vector<T> th = seq<T>(0.2, 2.5 / 4, Sim.nInitialConditions);
    //std::vector<size_t> nCollisions(Sim.nInitialConditions, 0);
    //
    std::vector<size_t> Indices = getIndicesOfVectorValue(Sim.TimesToGet, Sim.Time);
    std::cout << "Collecting displacements at times t = " << index(Sim.Time, Indices) << " with indices = " << Indices << std::endl;

    size_t M = Sim.nInitialConditions; //
    size_t count = 0;

    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < M; ++i)
    {
	bool ThreshCond = true;
        if (10 * count % Sim.nInitialConditions == 0 && count != 0 && ThreshCond)
            std::cout << "Simulated " << count << " initial conditions." << std::endl;

        // ON P.B.
        //Particle<T> P0 = Table.getIC_FixedLine(0, -1, Sim.TotalTime, Point<T>(-dx, (T)0), Point<T>(-dx, d));
        //P0.Vel = Point<T>(1,0);
        //
        //Particle<T> P0 = Table.getIC_FixedPosition(0, -1, Sim.TotalTime, Point<T>(0, dy + 0.0001));
        //P0.Vel = Point<T>(cos(th[i]), sin(th[i]));
        //std::cout << P0 << std::endl;
        
        // Standard distribution
        Particle<T> P0 = Table.getInitialState(0, -1, Sim.TotalTime);
        Particle<T> Pn = P0;

        Trajectory<T> Traj = Table.getTrajectory(Pn, Sim.Time);

        // Rejection sample
        std::vector<T> xDis = Traj.getDisplacements_x();
        for (size_t j = 0; j < Indices.size(); ++j) {
            T v = abs(xDis[Indices[j]]) / Sim.Time[Indices[j]]; // Make sure v = .5 at all times

	    //if (v < (Vbal - 0.005) && v > (Vbal + 0.005)) { // Only one instance of dx/t < .2 removes the sample
	    if (v < Vbal) {
		ThreshCond = false;
		break;
	    }
	}

	if (ThreshCond) {
	    for (size_t j = 0; j < Indices.size(); ++j) 
		Writer_DX << xDis[Indices[j]] << " ";
	    Writer_DX << std::endl;

	    // Write initial conditions
            //WriteVector<T>(Sim.DataDestination + "IC.dat", {P0.Pos.x, P0.Pos.y, P0.Vel.x, P0.Vel.y});
	    count++;
        } else {
	    M++;
	}
    }
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";
}



template <typename T>
inline void FirstPassageTimes(Billiard<T> &Table, Simulation<T> &Sim, const T& alpha)
{
    std::vector<T> MinFPT(Sim.L.size(), INFINITY); // Smallest first passage times
    Writer<T> _FPT(Sim.DataDestination + "FPT.dat", true);
    Writer<int> _FPC(Sim.DataDestination + "FPC.dat", true);
    Writer<T> _IC(Sim.DataDestination + "IC.dat", true);

    T PI = boost::math::constants::pi<T>();
    T dy = 1/(2 * tan(alpha/(T)2));
    T eps = (PI/(T)180);
    T beta = (PI - alpha)/(T)2;
    T beta_n = -beta + PI/(T)2;
    T d = Table.Obstacles[0]->getVertices()[5].y;
    //std::vector<T> y0 = seq<T>(1e-10, d, Sim.nInitialConditions); 

    T th0 = 16.3425218771 * PI/180;
    //T th0 = 15.42 * PI/180;
    //T th0 = 18.601 * PI/180;

    //std::vector<T> theta = seq<T>(-beta+(PI/(T)2)+eps, -beta+(PI/(T)2)-eps, Sim.nInitialConditions); // Bouncing ball
    //std::vector<T> theta = seq<T>(0, eps, Sim.nInitialConditions); // Propagating
    //std::vector<T> theta = seq<T>(th0-eps, th0, Sim.nInitialConditions); // Propagating
    std::vector<T> theta = seq<T>(th0-eps, th0, Sim.nInitialConditions); // Propagating
    //std::vector<T> theta = seq<T>(0, 2*PI, Sim.nInitialConditions); 

    // L = 729 + 1/2, y0 = d/2, x0=-dx + 1e-4, 
    //std::vector<T> theta = seq<T>(0.00015, 0.00017, Sim.nInitialConditions); // Propagating
    //std::vector<T> theta = seq<T>(0.0001566, 0.0001585, Sim.nInitialConditions); // Propagating

    Random<T> h0, phi0;

    std::cout << "d = " << d << ", dy = " << dy << ", dcos(alpha)+dy = " << d*cos(alpha)+dy << std::endl;
    //std::vector<T> theta = { beta + eps };
    std::vector<int> nUpdates(Sim.nInitialConditions, 0);

    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < Sim.nInitialConditions; ++i)
    {
        if (10 * i % Sim.nInitialConditions == 0 && i != 0)
            std::cout << "Simulated " << i << " initial conditions." << std::endl;

        Particle<T> IC = Table.getInitialState(0, -1, Sim.TotalTime);
	//T phi = phi0.getUniformRandom(beta_n - eps, beta_n + eps);
	//T phi = phi0.getUniformRandom(-eps, eps);
	//T phi = phi0.getUniformRandom(0, 2*PI);
	//IC.Vel = Point<T>(cos(phi), sin(phi));
	//IC.Vel = Point<T>(cos(theta[i]), sin(theta[i]));
	IC.Pos.x = (-1/(T)2) + 0.0000001;
	//IC.Pos.y = y0[i];
	IC.Pos.y = dy+0.1;
	
	// Trajectory example
	//IC.Pos.y = 0.3728098694;
	IC.Vel = Point<T>(cos(0), sin(0));
	//IC.Vel = Point<T>(cos(th0), sin(th0));
        Particle<T> IC1 = IC;


	//std::vector<T> ics = { IC.Pos.x, IC.Pos.y, phi };
	std::vector<T> ics = { IC.Pos.x, IC.Pos.y, theta[i] };
	_IC.WriteRowVector(ics);

	// Gets the true escape time (continuous)
	std::pair<std::vector<T>, std::vector<int>> data = Table.getEscapeTimes(IC, Sim.L);
	std::vector<T> EscapeTimes = data.first;
	_FPT.WriteRowVector(EscapeTimes);
	_FPC.WriteRowVector(data.second);

        //for (size_t k = 0; k < EscapeTimes.size(); ++k)
        //    if (EscapeTimes[k] < MinFPT[k])
        //        MinFPT[k] = EscapeTimes[k];

	// IC.nUpdates = num collisions - 2
	//Trajectory<T> dTraj = Table.getTrajectory(IC1, 100);
	//std::vector<T> x = dTraj.getLiftedPositions_x();
	//std::vector<T> y = dTraj.getLiftedPositions_y();
	//std::vector<T> fps = dTraj.getLiftedFreePaths();
	//std::vector<T> t(fps.size(), 0);
	//std::partial_sum(fps.begin(), fps.end(), t.begin());
	//WriteColVector<T>(Sim.DataDestination + "t.dat", t);
	//WriteColVector<T>(Sim.DataDestination + "x.dat", x);
	//WriteColVector<T>(Sim.DataDestination + "y.dat", y);
    }

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    std::cout << "\033[1;32mSimulation run time: " + 
	std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + 
	std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + 
	std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" << std::endl;
}

template <typename T>
inline void ScatterCell(Billiard<T> &Table, Simulation<T> &Sim, const T& alpha)
{
    std::vector<T> MinFPT(Sim.L.size(), INFINITY); // Smallest first passage times
    Writer<T> _FPT(Sim.DataDestination + "FPT.dat", true);
    Writer<int> _FPC(Sim.DataDestination + "FPC.dat", true);
    Writer<T> _IC(Sim.DataDestination + "IC.dat", true);

    T PI = boost::math::constants::pi<T>();
    T dy = 1/(2 * tan(alpha/(T)2));
    T eps = (PI/(T)180);
    T beta = (PI - alpha)/(T)2;
    T beta_n = -beta + PI/(T)2;
    T d = Table.Obstacles[0]->getVertices()[5].y;

    //std::vector<T> y0 = seq<T>(1e-10, d, Sim.nInitialConditions); 
    //std::vector<T> theta = seq<T>(th0-eps, th0, Sim.nInitialConditions); // Propagating

    Random<T> h0, phi0;
    std::cout << "d = " << d << ", dy = " << dy << ", dcos(alpha)+dy = " << d*cos(alpha)+dy << std::endl;

    std::vector<int> nUpdates(Sim.nInitialConditions, 0);
    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < Sim.nInitialConditions; ++i)
    {
        if (10 * i % Sim.nInitialConditions == 0 && i != 0)
            std::cout << "Simulated " << i << " initial conditions." << std::endl;

        Particle<T> IC = Table.getInitialState(0, -1, Sim.TotalTime);
	IC.Pos.x = (-1/(T)2) + 0.0000001;
	IC.Pos.y = h0.getUniformRandom(0, d);

	// Trajectory example
	IC.Vel = Point<T>(cos(0), sin(0));
        Particle<T> IC1 = IC;
	std::vector<T> ics = { IC.Pos.x, IC.Pos.y, 0 };
	_IC.WriteRowVector(ics);

	if (IC.Pos.y > dy) { // Collisionless trajectory with speed 1
	    std::vector<T> EscapeTimes = {1};
	    _FPT.WriteRowVector(EscapeTimes);
	} else if (IC.Pos.y > 0 && IC.Pos.y < d*cos(alpha) + dy) {
	    // Gets the true escape time (continuous)
	    std::pair<std::vector<T>, std::vector<int>> data = Table.getEscapeTimes(IC, Sim.L);
	    std::vector<T> EscapeTimes = data.first;
	    _FPT.WriteRowVector(EscapeTimes);
	} else {
	    std::vector<T> EscapeTimes = {0};
	    _FPT.WriteRowVector(EscapeTimes);
	}

	// IC.nUpdates = num collisions - 2
	//Trajectory<T> dTraj = Table.getTrajectory(IC1, 5);
	//std::vector<T> x = dTraj.getLiftedPositions_x();
	//std::vector<T> y = dTraj.getLiftedPositions_y();
	//std::vector<T> fps = dTraj.getLiftedFreePaths();
	//std::vector<T> t(fps.size(), 0);
	//std::partial_sum(fps.begin(), fps.end(), t.begin());
	//WriteVector<T>(Sim.DataDestination + "t.dat", t);
	//WriteVector<T>(Sim.DataDestination + "x.dat", x);
	//WriteVector<T>(Sim.DataDestination + "y.dat", y);
    }

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    std::cout << "\033[1;32mSimulation run time: " + 
	std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + 
	std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + 
	std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" << std::endl;
}



template <typename T>
inline void SFPT(Billiard<T> &Table, Simulation<T> &Sim)
{
    std::vector<T> MinFPT(Sim.L.size(), INFINITY); // Smallest first passage times

    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < Sim.nInitialConditions; ++i)
    {
        if (10 * i % Sim.nInitialConditions == 0 && i != 0)
            std::cout << "Simulated " << i << " initial conditions." << std::endl;

        Particle<T> IC = Table.getInitialState(0, -1, Sim.TotalTime);
        std::vector<T> EscapeTimes = Table.getEscapeTimes(IC, Sim.L);

        for (size_t k = 0; k < EscapeTimes.size(); ++k)
            if (EscapeTimes[k] < MinFPT[k])
                MinFPT[k] = EscapeTimes[k];

        //IC.WriteState(Sim.DataDestination + "InitialConditions.dat");
    }
    WriteVector(Sim.DataDestination + "SFPT.dat", MinFPT);

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";
}

// @brief Trajectories are conditionally generated, requiring that |x(t)| > Sim.L.back()
// @Write -Gets the sums of x position (w.r.t. time) of finite and infinite continuous trajectories
template <typename T>
inline void VelocityCorrelations(Billiard<T> &Table, Simulation<T> &Sim)
{

    Ensemble<T> VAC(Sim.Time.size());   // Velocity auto correlation
    Ensemble<T> VAC2D(Sim.Time.size()); // Velocity auto correlation
    Ensemble<T> IVAC(Sim.Time.size());  // Integrated velocity auto correlation

    Ensemble<T> TimeVAC(Sim.Time.size()); // VAC in time average
    Ensemble<T> BirkVAC(Sim.Time.size()); // VAC with incidence angles

    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < Sim.nInitialConditions; ++i)
    {
        if (10 * i % Sim.nInitialConditions == 0 && i != 0)
            std::cout << "Simulated " << i << " initial conditions." << std::endl;

        // For finite trajectories we need discrete trajectories
        Particle<T> FCT = Table.getInitialState(0, -1, Sim.TotalTime);
        Trajectory<T> Traj = Table.getTrajectory(FCT, Sim.Time);

        VAC.Push(Traj.getVAC_x());
        VAC2D.Push(Traj.getVAC_xy());
        IVAC.Push(Traj.getIVAC_x());

        TimeVAC.Push(CumSum<T>(Traj.getVAC_x()));
        BirkVAC.Push(Traj.getBirkVAC_x());
    }

    VAC.Write(Sim.DataDestination + "VAC.dat");
    VAC2D.Write(Sim.DataDestination + "VAC2D.dat");
    IVAC.Write(Sim.DataDestination + "IVAC.dat");
    TimeVAC.Write(Sim.DataDestination + "TimeVAC.dat");
    BirkVAC.Write(Sim.DataDestination + "BirkVAC.dat");

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";
    Sim.Print(RunTime);
    Sim.Write();
}

// @brief Trajectories are conditionally generated, requiring that |x(t)| > Sim.L.back()
// @Write - All the usual metrics of a finite system
template <typename T>
inline void EscapeToAbsorbingBoundary(Billiard<T> &Table, Simulation<T> &Sim)
{
    Ensemble<T> MeanEscapeTime(Sim.L.size());
    std::vector<Ensemble<T>> FiniteMoments(Sim.nMoments, Ensemble<T>(Sim.Time.size()));
    std::vector<Ensemble<T>> InfiniteMoments(Sim.nMoments, Ensemble<T>(Sim.Time.size()));

    std::vector<T> InfiniteTimes = getExponentialTimes<T>(100000, Sim.Time.size());

    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < Sim.nInitialConditions; ++i)
    {
        if (10 * i % Sim.nInitialConditions == 0 && i != 0)
            std::cout << "Simulated " << i << " initial conditions." << std::endl;

        // For finite trajectories we need discrete trajectories
        Trajectory<T> FiniteTraj, ContinTraj;

        Particle<T> IC = Table.getInitialState(0, -1, Sim.TotalTime);
        Particle<T> PC = IC;

        FiniteTraj = Table.getFiniteTrajectory(IC, Sim.L.back());
        ContinTraj = Table.getTrajectory(PC, InfiniteTimes);

        if (Table.Failed)
        {
            Table.ErrorReport(IC, i);
            i--; // Want M initial conditions
        }

        EscapeData<T> ED(Sim.L, FiniteTraj);
        ED.setEscapeData(Sim.L);
        MeanEscapeTime.Push(ED.EscapeTimes);

        // Finite moments
        Trajectory<T> FiniteContinuousTraj = FiniteTraj.MakeContinuous(Sim.Time); // Defaults to the last L
        if (i < 10)
            FiniteContinuousTraj.Write(Sim.DataDestination + "Traj_" + std::to_string(i) + ".dat");

        //if (ED.EscapeTimes.back() > (T)(0))
        FillMoments(FiniteMoments, FiniteContinuousTraj, Sim.MomentStep, Sim.nMoments);
        FillMoments(InfiniteMoments, ContinTraj, Sim.MomentStep, Sim.nMoments);

        WriteVector(Sim.DataDestination + "FP_Time.dat", ED.EscapeTimes);
        WriteVector(Sim.DataDestination + "FP_Abscissa.dat", ED.EscapePosition_x);
        WriteVector(Sim.DataDestination + "FP_Ordinate.dat", ED.EscapePosition_y);
        WriteVector(Sim.DataDestination + "FP_Direction.dat", ED.EscapeDirection);
        FiniteTraj.MakeParticle(0).WriteState(Sim.DataDestination + "InitialConditions.dat");
    }
    WriteVector(Sim.DataDestination + "MFPT.dat", MeanEscapeTime.PullMeans());
    for (size_t i = 0; i < FiniteMoments.size(); ++i)
    {
        FiniteMoments[i].Write(Sim.DataDestination + "FiniteMoment_" + std::to_string(i) + ".dat");
        InfiniteMoments[i].Write(Sim.DataDestination + "InfiniteMoment_" + std::to_string(i) + ".dat");
    }

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";
    Sim.Print(RunTime);
    Sim.Write();
}


// @brief Writes displacements at times t[Indices] 
template <typename T>
inline void getFront(Billiard<T> &Table, Simulation<T> &Sim, const T& vbal)
{
    Writer<T> _DX(Sim.DataDestination + "DX.dat", true); // Displacements
    Writer<T> _FP(Sim.DataDestination + "FP.dat", true); // Free paths
    Writer<T> _IC(Sim.DataDestination + "IC.dat", true); // Initial conditions

    // Geometry 
    //T d  = Table.Obstacles[0]->getVertices()[5].y;
    //T dy = Table.Obstacles[0]->getVertices()[1].y;
    //T dx = Table.Obstacles[0]->getVertices()[0].x;

    T eps = 0.075;
    T Vbal = vbal;

    // Some angles
    //std::vector<T> th = seq<T>(0.2, 2.5 / 4, Sim.nInitialConditions);
    //std::vector<size_t> nCollisions(Sim.nInitialConditions, 0);
   
    std::vector<size_t> Indices = getIndicesOfVectorValue(Sim.TimesToGet, Sim.Time);
    std::cout << "Collecting displacements at times t = " << index(Sim.Time, Indices) << " with indices = " << Indices << std::endl;

    size_t M = Sim.nInitialConditions; //
    size_t count = 0;
    bool flag = false;

    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < M; ++i)
    {
	bool ThreshCond = true;
        if (10 * count % Sim.nInitialConditions == 0 && count != 0 && flag) {
            std::cout << "Simulated " << count << " initial conditions." << std::endl;
	    flag = false;
	}

        // ON P.B.
        //Particle<T> P0 = Table.getIC_FixedLine(0, -1, Sim.TotalTime, Point<T>(-dx, (T)0), Point<T>(-dx, d));
        //P0.Vel = Point<T>(1,0);
        //
        //Particle<T> P0 = Table.getIC_FixedPosition(0, -1, Sim.TotalTime, Point<T>(0, dy + 0.0001));
        //P0.Vel = Point<T>(cos(th[i]), sin(th[i]));
        //std::cout << P0 << std::endl;
        
        // Standard distribution
        Particle<T> P0 = Table.getInitialState(0, -1, Sim.TotalTime);
        Particle<T> Pn = P0;

        Trajectory<T> Traj = Table.getTrajectory(Pn, Sim.Time);
        std::vector<T> xDis = Traj.getDisplacements_x();

        // Rejection sample from n.h. of ballistic velocity @Vbal +/- @eps
        for (size_t j = 0; j < Indices.size(); ++j) {
            T v = abs(xDis[Indices[j]]) / Sim.Time[Indices[j]]; 
	    if (v < Vbal - eps || v > Vbal + eps) { // Only one instance of dx/t < .2 removes the sample
		ThreshCond = false;
		break;
	    }
	}

	if (ThreshCond) {
	
	    // Write displacements
	    _DX.WriteRowVector(xDis, Indices);

	    // Write initial conditions
            //_IC.WriteRowVector({P0.Pos.x, P0.Pos.y, P0.Vel.x, P0.Vel.y});


	    // Need to get collisions until some maximum time is reached and then determing

	    // 
	    Particle<T> Pc = P0, Pt = P0;
            Trajectory<T> Traj1 = Table.getTrajectory(Pc, Sim.nCollisions);
            //Trajectory<T> Traj1 = Table.getTrajectory(Pc, Sim.nCollisions, Sim.Time[Indices[0]], Sim.Time[Indices.back()-1]);
	    
	    std::vector<T> LFP = Traj1.getLiftedFreePaths();
	    _FP.WriteRowVector(LFP);

	    // Save position series for plotting trajectories
	    ///Trajectory<T> Traj2 = Table.getTrajectory(Pt, Sim.Time);
            //Traj1.WritePositions(Sim.DataDestination + "ColTraj-" + std::to_string(i) + ".dat");
            //Traj2.WritePositions(Sim.DataDestination + "ConTraj-" + std::to_string(i) + ".dat");

	    flag = true;
	    count++;
        } else {
	    M++;
	}
    }

    //VB_EnseAvg.Write(Sim.DataDestination + "eMeanDXT.dat");
    //WriteVector(Sim.DataDestination + "tMeanDXT.dat", VB_TimeAvg);

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";
    std::cout << RunTime << std::endl;
}




#endif // MISC_INCLUDED
