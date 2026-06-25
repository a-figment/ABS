#include <iostream>
#include <algorithm>
#include <chrono>
#include <thread>
#include "Billiard.h"
#include "Trajectory.h"
#include "SimpleTrajectory.h"
#include "OnlineStatistics.h"
#include "Misc.h"
#include "Helper.h"
#include "Ensemble.h"
#include "Simulation.h"
#include "DataKeys.h"
#include "Writer.h"
#include "EscapeData.h"
#include "CoarseGrain.h"
#include "Beam.h"

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


template <typename T>
int getFlowDirection(const T& x, const Point<T>& Vel, const Point<long int>& cid, const T& dy, const T& dx) {

    Point<T> LHSVector = Point<T>(std::abs(dx), dy);
    Point<T> RHSVector = Point<T>(std::abs(dx), -dy);
    // Vector on each half of the elementary cell parallel to walls moving along + x-axis
    T FlowDir = 0;
    if (x < cid.x) // > 0 in either case => moving forward
	FlowDir = Dot<T>(LHSVector, Vel);
    else 
	FlowDir = Dot<T>(RHSVector, Vel);
    return (int)(FlowDir / std::abs(FlowDir));
}

template <typename T>
std::vector<Particle<T>> generateICs_thresh(Billiard<T>& Table, 
	const std::vector<T>& TIME, const int& M, const int& N, 
	const Point<T>& minstep, const T& mid, const T& eps, 
	const std::string& DataDir) {


    // -------------------------------------------------------
    T dx = Table.Obstacles[0]->getVertices()[0].x;
    T dy = Table.Obstacles[0]->getVertices()[1].y;
    Point<T> diam = Table.getDiameter();
    // -------------------------------------------------------
    
    std::vector<Particle<T>> particles;
    particles.reserve(M);

    // Record steps
    Writer<T> w_JUMP(DataDir + "jump.dat", true);
    size_t nSteps = 0;
    size_t maxSteps = 1000000; // Jumps to write

    // Main loop
    while (particles.size() <= (size_t)M) {

	Particle<T> P = Table.getInitialState(0, -1, TIME.back());
	Particle<T> P0 = P; // Doesn't get updated

	// Jump time/disp. and number of dir changes
	std::vector<Point<T>> steps;
	steps.reserve(N);
	std::vector<size_t> CC;
	CC.reserve(N);

	T x_init = P.LPos(diam).x;
	int start_flow = getFlowDirection(x_init, P.Vel, P.CellID, dy, dx); // initial flow direction
	int dir = start_flow;
	size_t n = 0;

	bool in_range = false;
	while (n < N) { // ALL steps along the traj are not in the range
	//while (P.TotalTime < TIME.back()) { // ALL steps along the traj are not in the range
	    T lt = 0, lx = 0, lx_init = x_init;
	    int steps_in_flow = 0; // Num steps along direction of flow
				   
	    while (start_flow == dir) {
		x_init = P.LPos(diam).x; // Position 
		T tt = P.TotalTime;
		Table.Update(P);
		lt += (P.TotalTime - tt); 
		start_flow = getFlowDirection(P.LPos(diam).x, P.Vel, P.CellID, dy, dx);
		steps_in_flow++;
		n++;
	    }
	    lx = P.LPos(diam).x - lx_init;
	    x_init = P.LPos(diam).x;

	    T v = std::abs(lx) / lt;
	    // Record jumps in the range or 
	    //if (nSteps < maxSteps) {
	    //    if (lt > minstep.x && std::abs(lx) > minstep.y) { 
	    //        steps.push_back(Point<T>(lt, lx));
	    //        nSteps++;
	    //    }
	    //}

	    // Reject trajectories which have a step velocity in a particular range
	    if (v < mid + eps && v > mid - eps) { // Inside interval
		if (lt > minstep.x && std::abs(lx) > minstep.y) { // 
		    in_range = true;
		    break; // reject
		}
	    } 
	    else {
	        if (nSteps < maxSteps) {
	            if (lt > minstep.x && std::abs(lx) > minstep.y) { 
	        	steps.push_back(Point<T>(lt, lx));
	        	nSteps++;
	            }
	        }
	    }
	    dir = -dir;
	}

	if (!in_range)
	    particles.push_back(P0);

	std::vector<T> lx = getCoordFromVectorPoint(steps, "x");
    	std::vector<T> ly = getCoordFromVectorPoint(steps, "y");
    	w_JUMP.WriteColVectors(lx, ly);
    }
    std::cout << "Generated " << particles.size() << " initial conditions." << std::endl;
    return particles;
}


template <typename T>
std::vector<Particle<T>> generateICs(Billiard<T>& Table, 
	const std::vector<T>& TIME, const int& M, 
	const Point<T>& minstep, const T& mid, const T& eps, 
	const std::string& DataDir) {


    // -------------------------------------------------------
    T dx = Table.Obstacles[0]->getVertices()[0].x;
    T dy = Table.Obstacles[0]->getVertices()[1].y;
    Point<T> diam = Table.getDiameter();
    // -------------------------------------------------------
    
    std::vector<Particle<T>> particles;
    particles.reserve(M);

    // Record steps
    Writer<T> w_JUMP(DataDir + "jump.dat", true);
    std::vector<Point<T>> steps;
    steps.reserve(M);
    size_t nSteps = 0;
    size_t maxSteps = 1000; // Jumps to write

    // Main loop
    while (particles.size() <= (size_t)M) {
	if (10 * particles.size() % (size_t)M == 0 && particles.size() != 0)
	    std::cout << "Generated " << particles.size() << " initial conditions." << std::endl;

	Particle<T> P = Table.getInitialState(0, -1, TIME.back());
	Particle<T> P0 = P; // Doesn't get updated

	std::vector<Point<T>> steps;
	std::vector<size_t> CC;
	T x_init = P.LPos(diam).x;
	int start_flow = getFlowDirection(x_init, P.Vel, P.CellID, dy, dx); // initial flow direction
	int dir = start_flow;
	size_t n = 0;

	bool in_range = false;
	while (P.TotalTime < TIME.back()) { // ALL steps along the traj are not in the range
	    T lt = 0, lx = 0, lx_init = x_init;
	    int steps_in_flow = 0; // Num steps along direction of flow
				   
	    while (start_flow == dir) {
		x_init = P.LPos(diam).x; // Position 
		T tt = P.TotalTime;
		Table.Update(P);
		lt += (P.TotalTime - tt); 
		start_flow = getFlowDirection(P.LPos(diam).x, P.Vel, P.CellID, dy, dx);
		steps_in_flow++;
		n++;
	    }
	    lx = P.LPos(diam).x - lx_init;
	    x_init = P.LPos(diam).x;

	    T v = std::abs(lx) / lt;
	    if (nSteps < maxSteps) {
		if (lt > minstep.x && std::abs(lx) > minstep.y) { 
		    if (!((v < mid + eps) && (v > mid - eps))) {
			steps.push_back(Point<T>(lt, lx));
			nSteps++;
		    }
		}
	    }

	    // Reject trajectories which have a step velocity in a particular range
	    if (v < mid + eps && v > mid - eps) {
		in_range = true;
		break;
	    } 
	    dir = -dir;
	}

	if (!in_range)
	    particles.push_back(P0);

	std::vector<T> lx = getCoordFromVectorPoint(steps, "x");
    	std::vector<T> ly = getCoordFromVectorPoint(steps, "y");
    	w_JUMP.WriteColVectors(lx, ly);
    }
    std::cout << "Generated " << particles.size() << " initial conditions." << std::endl;
    return particles;
}

// @brief Writes moments of displacement
template <typename T>
inline void getMoments_t(Billiard<T> &Table, Simulation<T> &Sim, const std::vector<Particle<T>>& Ps)
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

	Particle<T> P0;
        if (Ps.size() == (size_t)Sim.nInitialConditions+1) {
	    P0 = Ps[i];
        } else {
	    P0 = Table.getInitialState(0, -1, Sim.TotalTime);
        }

//	Particle<T> P0 = Table.getInitialState(0, -1, Sim.TotalTime);
//	T th = atan2(P0.Vel.y, P0.Vel.x);
//	while (in_nh<T>(th, 0, eps) || 
//	       in_nh<T>(th, PI, eps) || 
//	       in_nh<T>(th, -alpha, eps) || 
//	       in_nh<T>(th, PI - alpha, eps)) {
//	    P0 = Table.getInitialState(0, -1, Sim.TotalTime);
//	    th = atan2(P0.Vel.y, P0.Vel.x);
//	}
//        Ps.push_back(P0);
	//Particle<T> P0 = Table.getInitialState(0, -1, Sim.TotalTime);
	//P0.Vel = Vels[i];

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
	Moments[i].Write(Sim.DataDestination + "Moment_" + std::to_string(i) + ".dat");
}

// --------------------------------------------------------------------
// -------------------------------- MAIN ------------------------------
// --------------------------------------------------------------------
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
    DataKeys DK({"p", "a", "d", "i", "c"}, {1000, 40, 40, 6, 3}, {"PARTITIONS", "ALPHA_R = PI/2", "PORE WIDTH", "STEP ID", "STEP COMBINATION"});
    std::string Key = DK.MakeKey(KeyIndices);
    DK.Print();
    // ======================= BOILER PLATE END ==========================

    // ======================= BILLIARD DEFINITION START ==========================
    std::vector<std::unique_ptr<Obstacle<float_>>> Obstacles(1);
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
    std::vector<float_> ratio = seq<float_>(1/(float_)4, 1, 25); // r = 1/4, 1/2, 3/4, 1 -> 0, 8, 16, 25
    std::vector<float_> porewidths(ratio.size(), 0); //= { - 1/((float_)2), (float_)1 };
    for (size_t i = 0; i < porewidths.size(); ++i) {
        porewidths[i] = ratio[i] * dy;
    }
    Polygon<float_> Channel = getParallelChannel2<float_>(dy, (1 / (float_)2), porewidths[KeyIndices[2]]);
    //Polygon<float_> Channel = getParallelChannel2<float_>(3/(float_)10, (float_)1, 1/(float_)10);
    //Polygon<float_> Channel = getParallelChannel2<float_>((float_)10, (float_)1, 51/(float_)100);
    Obstacles[0] = std::make_unique<Polygon<float_>>(Channel);
    Billiard<float_> Table(Obstacles);
    // ======================= BILLIARD DEFINITION END ==========================



    // =============================== SIMULATION START ===============================
    std::string DataDir = "Data/" + Key;
    //std::string DataDir = "/fred/oz055/jorchard_data/ABS/S-MID/" + Key;
    
    const float_ TotalTime = 1e5;
    std::vector<float_> TIME = seq<float_>(0, TotalTime, TotalTime/1000 + 1);

    Simulation<float_> Sim {
            .Time = TIME,
            .L = {0},
            .TimesToGet = {1e4-1, 1e5-1, 1e6-1},
            .Symbols = {""},
            .DataDestination = DataDir,
            .RunTime = ".",
            .MomentStep = (float_)1,
            .TotalTime = TIME.back(),
            .nInitialConditions = 100000,
            .nCollisions = 1000000, // always reach max time before max coll
            .nMoments = 2};

    // Neighborhoods for rational and irrational special velocities
    // PI/3, PI/2, 2PI/3
    float_ aeps = 25 / (float_)1000;
    std::vector<float_> mids_rat = { 19/(float_)100, 1/(float_)2, 6/(float_)10 };
    std::vector<float_> eps_rat  = { aeps, aeps, aeps };
    std::vector<float_> mids_irr = { 31/(float_)100, 1/(float_)2, 6/(float_)10 };
    std::vector<float_> eps_irr  = { aeps, aeps, aeps };
    
    float_ mid = 0, eps = 0;
    switch (KeyIndices[1]) {
	case 6: 
	    mid = mids_rat[0];
	    eps = eps_rat[0];
	    break;
	case 7:
	    mid = mids_rat[1];
	    eps = eps_rat[1];
	    break;
	case 8:
	    mid = mids_rat[2];
	    eps = eps_rat[2];
	    break;
	case 21:
	    mid = mids_irr[0];
	    eps = eps_irr[0];
	    break;
	case 22:
	    mid = mids_irr[1];
	    eps = eps_irr[1];
	    break;
	case 23:
	    mid = mids_irr[2];
	    eps = eps_irr[2];
	    break;
    }

    std::vector<float_> mintimes = { -1, (float_)1, (float_)10, (float_)100, (float_)1000, (float_)10000  };
    std::vector<float_> mindists = { -1, (float_)1, (float_)10, (float_)100, (float_)1000, (float_)10000  };
    size_t step_id = KeyIndices[3];
    std::vector<Point<float_>> minsteps = { Point<float_>(mintimes[step_id], -1), Point<float_>(-1, mindists[step_id]), Point<float_>(mintimes[step_id], mindists[step_id])};

    std::cout << "Jump bound = " << minsteps[KeyIndices[4]] << ". Neighbourhood = [" << mid << " - " << eps << ", " << mid << " + " << eps << "]. Aperture angle = " << alpha[KeyIndices[1]] * 180/PI << std::endl;

    std::vector<Particle<float_>> Ps;
    if (KeyIndices[4] <= 2) {
	std::cout << "Generating modified initial distribution." << std::endl;
	//Ps = generateICs<float_>(Table, TIME, Sim.nInitialConditions, minsteps[KeyIndices[4]], mid, eps, DataDir);
	Ps = generateICs_thresh<float_>(Table, TIME, Sim.nInitialConditions, Sim.nCollisions, minsteps[KeyIndices[4]], mid, eps, DataDir);
	WriteParticles<float_>(DataDir + "P0.dat", Ps);
    }
    //getMoments_t<float_>(Table, Sim, Ps);


    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";
    std::cout << RunTime << std::endl;
    return 0;
}

