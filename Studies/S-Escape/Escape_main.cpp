#include <iostream>
#include <algorithm>
#include <chrono>
#include <functional>
#include <thread>
#include "Billiard.h"
#include "Trajectory.h"
#include "OnlineStatistics.h"
#include "Misc.h"
#include "Helper.h"
#include "Simulation.h"
#include "DataKeys.h"
#include "Writer.h"
#include "Logger.h"
//#include "Scatter.h"
//#include "ScatterV2.h"
#include "ScatterV3.h"
#include "Ensemble.h"
// THIS TEMPLATE INTENDS TO RESET THE STANDARD FOR 'STUDIES'. THE GOAL IS THE AVOID
// USING THE MISC.H FUNCTIONS TO RUN BILLIARD SIMULATIONS UNLESS IT IS COMPLETELY GENERIC.
//
// EXAMPLE USAGE OF MISC.H FUNCTIONS: 
//		1. COLLECTING MOMENTS OF DISPLACEMENT
//		2. CORRELATION FUNCTIONS (POSITION, VELOCITY, INTEGRATED VELOCITY ETC.)
//		3. 

// < single character prefix, max number of instantiations >
typedef float_ldub float_;
//typedef float_50 float_;
const int nDigits = std::numeric_limits<float_>::digits10;

template<typename T>
std::vector<int> intf(const std::vector<T>& x) {

    std::vector<size_t> gz, lz;
    gz.reserve(x.size());
    lz.reserve(x.size());
    for (size_t i = 0; i < x.size(); i++) {
	if (x[i] >= 0)
	    gz.push_back(i);
	else
	    lz.push_back(i);
    }

    std::vector<int> val(x.size());
    if (gz.size() == 0) { // All negative
	for (size_t i = 0; i < x.size(); ++i) 
	    val[i] = floor(x[i]);
    } else {		  // All positive
	for (size_t i = 0; i < gz.size(); ++i) 
	    val[gz[i]] = ceil(x[gz[i]]);
	for (size_t i = 0; i < lz.size(); ++i) 
	    val[lz[i]] = floor(x[lz[i]]);
    }
    return val;
}

template<typename T>
T getCrossingHeight(const std::vector<T>& x, const std::vector<T>& y) {
    int mid = 0;
    T yint = 0;
    if (floor(x[0]) != floor(x[1])) {
	mid = std::max(floor(x[0]), floor(x[1]));
	T m = (y[1] - y[0])/(x[1] - x[0]);
	yint = m*(mid-x[1]) + y[1];
	return yint;
    } else {
	return INFINITY;
    }
}

//template<typename T>
//T getCrossingLocation(const std::vector<T>& x, const std::vector<T>& y) {
//    int mid = 0;
//    T yint = 0;
//    if (floor(x[0]) != floor(x[1])) {
//	mid = std::max(floor(x[0]), floor(x[1]));
//	return mid;
//    } else {
//	return INFINITY;
//    }
//}
template<typename T>
T getCrossingLocation(const std::vector<T>& x) {
    int mid = 0;
    if (floor(x[0]) != floor(x[1])) {
	mid = std::max(floor(x[0]), floor(x[1]));
	return mid;
    } else {
	return INFINITY;
    }
}

template<typename T>
std::vector<T> getCrossingHeights(const std::vector<T>& x, const std::vector<T>& y) {

    std::vector<T> Heights;
    Heights.reserve(x.size());
    Heights.push_back(y[0]);
    for (size_t i = 0; i < x.size()-1; ++i) {
	std::vector<T> X = {x[i], x[i+1]};
	std::vector<T> Y = {y[i], y[i+1]};
	T F = getCrossingHeight(X, Y);
	if (F != INFINITY)
	    Heights.push_back(F);
    }
    return Heights;
}

template<typename T>
std::vector<int> getCrossingLocations(const std::vector<T>& x, const std::vector<T>& y) {

    std::vector<int> Location;
    Location.reserve(x.size());
    Location.push_back(0);
    for (size_t i = 0; i < x.size()-1; ++i) {
	std::vector<T> X = {x[i], x[i+1]};
	std::vector<T> Y = {y[i], y[i+1]};
	int F = getCrossingLocation(X);
	if (F != std::numeric_limits<int>::max())
	    Location.push_back(F);
    }
    return Location;
}


template<typename T>
void NoRecursion(const T& d, const T& eps, const size_t& M,  const std::vector<T>& ttg) {
    Random<T> UR;

    // WITHOUT RECURSION
    for (size_t k = 0; k < M; ++k) {

	T Angle = UR.getUniformRandom(0, 2*PI - eps);
	T Height = UR.getUniformRandom((T)0 + eps, d - eps);

	// Initial conditions
	T h0 = Height;
	T T0 = Angle;
	T tau = 0;
	int_ll x0 = 0;
	int bid = -1;
	bool sym = false;

	size_t tidx = 0;
	T TAU = tau;

	Scatterer<T> S_(d, h0, T0, x0, bid);
	while (TAU < ttg.back()) {
	    S_.Evolve_(h0, T0, x0, tau, bid, sym);
	    sym = false;
	    TAU += tau;
	    if (TAU >= ttg[tidx] && tidx < ttg.size()) {
		tidx++;
	    }
	}
    }
}

template<typename T>
std::vector<T> getHeights(const size_t& nToGet, const T& eps, const T& d) {
    Random<T> UR;
    std::vector<T> Heights(nToGet, -1);
    for (size_t i = 0; i < nToGet; ++i) 
	Heights[i] = UR.getUniformRandom(eps, d - eps);
    return Heights;
}

template<typename T>
std::vector<T> getAngles(const size_t& nToGet, const T& eps) {
    Random<T> UR;
    std::vector<T> Angles(nToGet, -1);
    for (size_t i = 0; i < nToGet; ++i)  {
	//Angles[i] = UR.getUniformRandom(0, 2*PI - eps);
	Angles[i] = UR.getUniformRandom(-PI/(float_)2 + eps, (3*PI/(float_)2) - eps);
	//Angles[i] = UR.getUniformRandom(0, (PI/(float_)4) - eps);
    }
    return Angles;
}

// Angles between (3pi/2, 2pi] U [0, pi/2)
// @nToGet must be an even number
template <typename T>
std::vector<T> getAnglesRHP(const size_t& nToGet, const T& eps) {
    std::cout << "ANGLES NOT FILLING [0, 2pi]" << std::endl;
    Random<T> UR1;
    Random<T> UR2;
    std::vector<T> Angles(nToGet, -1);
    for (size_t i = 0; i < nToGet/2; ++i) 
	Angles[i] = UR1.getUniformRandom((3*PI/(T)2) + eps, 2*PI - eps);
    for (size_t i = nToGet/2; i < nToGet; ++i) 
	Angles[i] = UR2.getUniformRandom(0, (PI/(T)2) - eps);
    return Angles;
}

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
    DataKeys DK({"p", "a", "r"}, {1000, 100, 100}, {"PARTITIONS", "ALPHA_R = PI/4, PI/3, ..., 3PI/4", "RATIO"});
    std::string Key = DK.MakeKey(KeyIndices);
    DK.Print();

    // Consts
    float_ PHI = boost::math::constants::phi<float_>();
    float_ PI = boost::math::constants::pi<float_>();
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
    //float_ alpha_ = 2*PI/(float_)3;
    //float_ dy = dx * cos(alpha_ / (float_)2) / sin(alpha_ / (float_)2);

    // PARALLEL CHANNEL
    //std::vector<float_> D = { 1/(float_)8, 2/(float_)5, 1/(float_)2, (float_)1 };
    //std::vector<float_> ratio = seq<float_>(1, abs(1/(cos(alpha_))), 50);
    //std::vector<float_> ratio = seq<float_>(1, abs(1/(cos(alpha[KeyIndices[1]]))), 100);
    std::vector<float_> ratio = seq<float_>(1, 10, 20);
    float_ d = ratio[KeyIndices[2]] * dy;

    //Polygon<float_> Channel({0, dx, 2*dx, 2*dx, dx, 0}, {0, dy, 0, d, d + dy, d}, {2, 5}, false);
    Polygon<float_> Channel({-dx, 0, dx, dx, 0, -dx}, {0, dy, 0, d, d + dy, d}, {2, 5}, false);
    std::cout << Channel << std::endl;

    Obstacles[0] = std::make_unique<Polygon<float_>>(Channel);
    Billiard<float_> Table(Obstacles);

    // =============================== SIMULATION START ===============================    
    std::string DataDir = "Data/" + Key;
    //std::string DataDir = "/fred/oz055/pathtohpc/ABS/Scatter/" + Key;

    float_ eps = 1e-12;
    size_t M = 100000;

    std::vector<float_> Time = seq<float_>(0, 5000, 1001);

    Simulation<float_> Sim {
        .Time = Time,
        .L = { 1+dx },
        //.L = { dx, 3 + dx, 4 + dx, 5 + dx },
        //.L = { 9 + dx, 27+dx, 81+dx },//, 27 + dx }, //81 + dx, 243 + dx, 729 + dx },
        .TimesToGet = {4000, 8000 },
        .Symbols = {""},
        .DataDestination = "Data/" + Key,
        .RunTime = ".",
        .MomentStep = 2,
        .TotalTime = Time.back(),
        .nInitialConditions = M,
        .nCollisions = 0, // always reach max time before max coll
        .nMoments = 1};
    getDisplacements_t<float_>(Table, Sim, {551, 751, 1000});
    //FirstPassageTimes<float_>(Table, Sim, alpha[KeyIndices[1]]);
    //ScatterCell<float_>(Table, Sim, alpha_);
    //ScatterCell<float_>(Table, Sim, alpha[KeyIndices[1]]);

    //Particle<float_> P0 = Table.getInitialState(0, -1, 1e7);
    //float_ he = Heights[0];
    //float_ th = Angles[0];
    //P0.Pos = Point<float_>(0, he);
    //P0.Vel = Point<float_>(cos(th), sin(th));
    //

    //Trajectory<float_> dTraj = Table.getTrajectory(P0, 100);
    //std::vector<float_> x = dTraj.getLiftedPositions_x();
    //std::vector<float_> y = dTraj.getLiftedPositions_y();
    //std::vector<float_> fps = dTraj.getLiftedFreePaths();
    //std::vector<float_> t(fps.size(), 0);
    //std::partial_sum(fps.begin(), fps.end(), t.begin());
    ////WriteColVector<float_>(DataDir + "t.dat", t);
    //WriteColVector<float_>(DataDir + "x.dat", x);
    //WriteColVector<float_>(DataDir + "y.dat", y);

    //std::vector<int> LocationFromTraj = getCrossingLocations(x, y);
    //std::vector<float_> HeightsFromTraj = getCrossingHeights(x, y);
    //WriteVector<float_>(DataDir + "RawH.dat", HeightsFromTraj);
    //WriteVector<int>(DataDir + "RawX.dat", LocationFromTraj);
    // ====================================================


    // COMPARING PERFORMANCE OF RECURSIVELY AND NON-RECURSIVELY DEFINED MAPS
    //size_t M= 10000; // Number of initial conditions
    //size_t N = 100; // Number of runs
    //std::vector<float_> ttg = seq<float_>(1, 1e4, 1000);
    //std::vector<int> rtime(N,0);
    //std::vector<int> nrtime(N,0);
    //for (size_t i = 0; i < N; ++i) {
    //    if (10 * i % N == 0 && i != 0)
    //        std::cout << "Simulated " << i << " runs." << std::endl;

    //    std::chrono::high_resolution_clock::time_point t0_recursion = std::chrono::high_resolution_clock::now();
    //    Recursion<float_>(d, eps, M, ttg);
    //    std::chrono::high_resolution_clock::time_point t1_recursion = std::chrono::high_resolution_clock::now();
    //    rtime[i] = std::chrono::duration_cast<std::chrono::milliseconds>(t1_recursion - t0_recursion).count();


    //    std::chrono::high_resolution_clock::time_point t0_normal = std::chrono::high_resolution_clock::now();
    //    NoRecursion<float_>(d, eps, M, ttg);
    //    std::chrono::high_resolution_clock::time_point t1_normal = std::chrono::high_resolution_clock::now();
    //    nrtime[i] = std::chrono::duration_cast<std::chrono::milliseconds>(t1_normal - t0_normal).count();
    //}
    //std::cout << rtime << nrtime << std::endl;
    //WriteVector<int>(DataDir + "RunTimes.dat", rtime);
    //WriteVector<int>(DataDir + "RunTimes.dat", nrtime);

//    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
//    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";
//    std::cout << RunTime << std::endl;
    return 0;
}


    //std::vector<float_> TIME = seq<float_>(1, 1e7+1, 10001);
    //std::vector<float_> ttg = TIME;
    //ttg.erase(ttg.begin());
    //Simulation<float_> Sim {
    //    .Time = TIME,
    //    .L = {0},
    //    .TimesToGet = ttg,
    //    .Symbols = {""},
    //    .DataDestination = "Data/" + Key,
    //    .RunTime = ".",
    //    .MomentStep = 0,
    //    .TotalTime = TIME.back(),
    //    .nInitialConditions = 2000,
    //    .nCollisions = 1, // always reach max time before max coll
    //    .nMoments = 0};
    //Sim.Write(); // Write the simulation object?
    

