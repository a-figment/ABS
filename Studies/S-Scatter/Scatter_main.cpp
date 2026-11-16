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
#include "ScatterV2.h"
//#include "ScatterV3.h"
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


//template<typename T>
//std::vector<int> intf(const std::vector<T>& x) {
//
//    std::vector<size_t> gz, lz;
//    gz.reserve(x.size());
//    lz.reserve(x.size());
//    for (size_t i = 0; i < x.size(); i++) {
//	if (x[i] >= 0)
//	    gz.push_back(i);
//	else
//	    lz.push_back(i);
//    }
//
//    std::vector<int> val(x.size());
//    if (gz.size() == 0) { // All negative
//	for (size_t i = 0; i < x.size(); ++i) 
//	    val[i] = ceil(x[i]);
//    } else {		  // All positive
//	for (size_t i = 0; i < gz.size(); ++i) 
//	    val[gz[i]] = floor(x[gz[i]]);
//	for (size_t i = 0; i < lz.size(); ++i) 
//	    val[lz[i]] = ceil(x[lz[i]]);
//    }
//
//
//    for (size_t i = 0; i < val.size(); ++i) {
//	if (val[i] == 0 && x[i] > 0)
//	    val[i] = +1;
//	if (val[i] == 0 && x[i] < 0)	
//	    val[i] = -1;
//    }
//
//    return val;
//
//}

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

//template<typename T>
//void Recursion(const T& d, const T& eps, const size_t& M,  const std::vector<T>& ttg) {
//    Random<T> UR;
//    for (size_t k = 0; k < M; ++k) {
//
//	float_ Angle = UR.getUniformRandom(0, 2*PI - eps);
//	float_ Height = UR.getUniformRandom((float_)0 + eps, d - eps);
//
//	// Initial conditions
//	float_ h0 = Height;
//	float_ T0 = Angle;
//	float_ tau = 0;
//	int_ll x0 = 0;
//	int bid = -1;
//	bool sym = false;
//
//	size_t tidx = 0;
//	float_ TAU = tau;
//
//	Scatterer<float_> S(d, h0, T0, x0, bid);
//	while (TAU < ttg.back()) {
//	    S.Evolve(h0, T0, x0, tau, bid, sym);
//	    sym = false;
//	    TAU += tau;
//	    if (TAU >= ttg[tidx] && tidx < ttg.size()) {
//		tidx++;
//	    }
//	}
//    }
//}

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

template<typename T>
std::vector<T> getAnglesRHS(const size_t& nToGet, const T& eps) {
    Random<T> UR;
    std::vector<T> Angles(nToGet, -1);
    for (size_t i = 0; i < nToGet; ++i)  {
	Angles[i] = UR.getUniformRandom(-PI/(float_)2 + eps, (PI/(float_)2) - eps);
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


// @brief Takes a uniformly random initial condition and checks if it is in the @branchToGet
//	  NOTE: This does not ensure that the regions are uniformly populated
//	  NOTE: There is no exit condition. If the loop does not terminate, it means there is a 
//		mistake in the region definition
// @branchToGet An integer between 0-18 
template <typename T>
Point<T> getPointInBranch(const int& branchToGet, const T& eps) {
    Random<T> URA;
    Random<T> URH;

    // One iteration per IC
    size_t n = 0;
    T d = 1/(T)2;
    T h0 = URH.getUniformRandom(eps, d - eps);
    //T T0 = URA.getUniformRandom(0, 2*PI - eps);
    T T0 = URA.getUniformRandom((-PI/(T)2) + eps, (3*PI/(T)2) - eps);
    T tau = 0;
    int_ll x0 = 0;
    int bid = -1;
    bool sym = false;

    T foundH0 = -1;
    T foundT0 = 0;
    do {
	h0 = URH.getUniformRandom(eps, d - eps);
	//T0 = URA.getUniformRandom((3*PI/(T)2) + eps, 2*PI);
	//T0 = URA.getUniformRandom(0, 2*PI - eps);
	T0 = URA.getUniformRandom((-PI/(T)2) + eps, (3*PI/(T)2) - eps);
	foundH0 = h0;
	foundT0 = T0;
	tau = 0;
	x0 = 0;
	bid = -1;
	sym = false;

	Scatterer<T> S(d, h0, T0, x0, bid);	// Evolves the map
	S.Evolve_(h0, T0, x0, tau, bid, sym);

	n++;
    } while(branchToGet != bid); 

    return Point<T>(foundH0, foundT0);
}

// @brief Returns vector of height and angles which belong to branch @ID
template <typename T>
std::pair<std::vector<T>, std::vector<T>> getPointsInBranch(const int& ID, const T& eps, const T& d, const size_t& nPoints) {

    std::cout << "Getting points in branch " << ID << std::endl;
    std::vector<T> Heights(nPoints,0);
    std::vector<T> Angles(nPoints,0);
    for (size_t j = 0; j < nPoints; ++j) {
	Point<T> P = getPointInBranch<T>(ID, eps);
	Heights[j] = P.x;
	Angles[j] = P.y;
    }
    WriteVector<T>("Data/Branch-" + std::to_string(ID) + ".dat", Heights);
    WriteVector<T>("Data/Branch-" + std::to_string(ID) + ".dat", Angles);
    return std::make_pair(Heights, Angles);
}


// @brief Evolves @M particles until they change direction and records the time/distance
//	  travelled. Only the first "jump" is considered before that trajectory is forgotten
// @minTime Is the smallest time to record
template <typename T>
void JumpTimes(const size_t& M, const T& minTime, const T& d, const std::string& DataDir) {
    
    T eps = 1e-12;
    Random<T> URH, URA;

    size_t nJumps = 0;

    Logger<T> LOG_ics(DataDir + "ICs.dat", true);

    std::vector<T> JumpTimes(M, -1);
    std::vector<T> JumpDistances(M, -1);
    while (nJumps < M) {
	//if (10 * nJumps % M == 0 && nJumps != 0)
	//    std::cout << "Simulated " << nJumps << " initial conditions." << std::endl;

	T h0 = URH.getUniformRandom(eps, d - eps);
	T T0 = URA.getUniformRandom(0, 2*PI - eps);

	T tau = 0;
	int_ll x0 = 0;
	int bid = -1;
	bool sym = false;

	size_t j = 0;
	T TAU = 0;

	Scatterer<T> S(d, h0, T0, x0, bid);	// Evolves the map
	int_ll xPrev = x0;
	do {
	    xPrev = x0;
	    TAU += tau;
	    S.Evolve_(h0, T0, x0, tau, bid, sym);
	    //X.push_back(x0);
	    sym = false;
	    //std::cout << Point<int_ll>(xPrev, x0) << " " << TAU << " " << TAU + tau << " " << nJumps << std::endl;
	} while (abs(xPrev) < abs(x0));

	//std::cout << X << std::endl;
	if (TAU > minTime) {
	    LOG_ics.Stream << h0 << " " << T0 << std::endl;
	    JumpTimes[nJumps] = TAU;
	    JumpDistances[nJumps] = xPrev;
	    nJumps++;
	}
    }

    WriteColVector<T>(DataDir + "JumpTimes.dat", JumpTimes);
    WriteColVector<T>(DataDir + "JumpDistances.dat", JumpDistances);
}

// @brief Writes trajectories data which is a particular speed
template <typename T>
void SearchAndWrite(const std::string& filename, const std::vector<T>& speeds, const std::vector<T>& lbs, const std::vector<T>& ubs, const T& h0, const T& t0, const T& d, const bool& WriteTrajectory) {

    assert(speeds.size() == 2);
    for (size_t i = 0; i < lbs.size(); ++i) {
	bool inRange = in_nh_ublb(abs(speeds[0]), lbs[i], ubs[i]) && in_nh_ublb(abs(speeds[1]), lbs[i], ubs[i]);
	if (inRange) { // Speed was inside range - now write trajectory data

	    Logger<T> w_VBLOG(filename + std::to_string(i) + ".dat", true);

	    if (WriteTrajectory) {
		// Trajectory details - Actually simulate another trajectory and extract data
		size_t N = 50000;
		std::vector<size_t> itg = seq<size_t>(0, N, N+1);
		ScatterTraj<T> STn = getTrajectory(h0, t0, d, itg);
		// Write initial state and all other states
		int branch = whichBranch<float_>(h0, t0, d);
		w_VBLOG.Stream << h0 << " " << t0 << " " << 0 << " " << 0 << " " << branch << std::endl;
		for (size_t i = 1; i < itg.size(); ++i) {
		    w_VBLOG.Stream << STn.Height[i] << " " << STn.Angle[i] << " " << 
			STn.Position[i] << " " << STn.Tau[i] << " " << STn.BranchID[i] << std::endl;
		}
	    } else {
		// Just write the initial conditions
		int branch = whichBranch<T>(h0, t0, d);
		std::cout << std::fixed << std::setprecision(5) << "writing! " << i << ": " << lbs[i] << " < (" << speeds[0] << ", " << speeds[1] << ") < " << ubs[i] << std::endl;
		w_VBLOG.Stream << h0 << " " << t0 << " " << branch << " " << speeds[0] << " " << speeds[1] << std::endl;
	    }
	}
    }
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
    DataKeys DK({"p", "a", "r"}, {1000, 10, 100}, {"PARTITIONS", "ALPHA_R = PI/4, PI/3, ..., 3PI/4", "RATIO"});
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

    // PARALLEL CHANNEL
    std::vector<float_> D = { 1/(float_)2, (float_)1 };
    float_ d = D[KeyIndices[2]];
    Polygon<float_> Channel({0, dx, 2*dx, 2*dx, dx, 0}, {0, dy, 0, d, d + dy, d}, {2, 5}, false);
    std::cout << Channel << std::endl;

    Obstacles[0] = std::make_unique<Polygon<float_>>(Channel);
    Billiard<float_> Table(Obstacles);

    // =============================== SIMULATION START ===============================    
    std::string DataDir = "Data/" + Key;
    //std::string DataDir = "/fred/oz055/pathtohpc/ABS/S-Scatter/" + Key;

    float_ eps = 1e-12;
    size_t M = 1;

    // x - 315 = num found in XXX time
    
    // =============================== INITIAL CONDITIONS START ===============================    
    //std::vector<float_> SpecialAngles = { 0, PI, 
    //    PI/(float_)4, PI + PI/(float_)4, 
    //    (PI/(float_)2) + PI/(float_)4, (3*PI/(float_)2) + PI/(float_)4,
    //    atan(2), PI + atan(2), (PI/(float_)2) + atan(2), (3*PI/(float_)2) + atan(2),
    //    atan(1/(float_)2), PI + atan(1/(float_)2), (PI/(float_)2) + atan(1/(float_)2), (3*PI/(float_)2) + atan(1/(float_)2) };

    //size_t M = 10000+SpecialAngles.size(); // Number of initial conditions
    //std::vector<float_> Angles = getAngles(M-SpecialAngles.size(), eps);
    //for (size_t i = 0; i < SpecialAngles.size(); ++i)
    //    Angles.insert(Angles.begin(), SpecialAngles[i]);
    //std::vector<float_> Heights = getHeights(M, eps, d);

    // Uniform random IC's (the DEFAULT)
    //std::vector<float_> Angles = getAngles(M, eps);
    //std::vector<float_> Heights = getHeights(M, eps, d);

    // A grid of IC's which cover D \times A
    //float_ eps1 = 0.0001;
    //std::vector<float_> Heights = seq<float_>(eps1, d-eps1, M);
    //std::vector<float_> Angles = seq<float_>((-PI/(float_)2) + eps1, (PI/(float_)2) - eps1, std::round(M));

    // ============================================================================
    // Write region data - Takes angles in RHS unit circle and updates the map ONCE
    //std::vector<float_> Angles = getAnglesRHS(M, eps);
    //Logger<float_> PointLog(DataDir + "SV.dat", true);
    //for (size_t i = 0; i < Angles.size(); ++i) {
    //    float_ h0 = Heights[i];
    //    float_ T0 = Angles[i];
    //    float_ tau = 0;
    //    int_ll x0 = 0;
    //    int bid = -1;
    //    bool sym = false;
    //    Scatterer<float_> S(d, h0, T0, x0, bid)	    // Evolves the map
    //    S.Update(h0, T0, x0, tau, bid, sym);	    // Accepts only angles in RHS circle

    //    PointLog.Stream << Heights[i] << " " << Angles[i] << " " << tau << " " << x0 << " " << bid << std::endl;
    //    //`ScatterVar<float_> SV = getBranchPoint(d, Heights[i], Angles[i]);
    //    //PointLog.Stream << SV.Height << " " << SV.Angle << " " << SV.Tau << " " << SV.Position << " " << SV.BranchID << std::endl;
    //}
    //return 0;
    // ============================================================================

    // User defined IC
    //std::vector<float_> Heights = seq<float_>(eps, d - eps, M); // { 7/(float_)100 };
    //std::vector<float_> Angles(M, (-PI/(float_)3)); // = { atan(1/(float_)2) };
    std::vector<float_> Heights = seq<float_>(eps, d - eps, M);
    std::vector<float_> Angles(M, -0.8); 
    //Heights = { 0.19943002912272900229, 0.19943002912272900229};
    //Angles = { -7.1309089151128518536e-02, 7.1309089151128518536e-02};
    float_ hhh =  0.19943002912;
    float_ aaa = -0.07130908915;
    Heights = { hhh, hhh };
    Angles = { aaa, -aaa };

    std::cout << Heights;
    std::cout << Angles;
    //std::vector<float_> Angles = seq<float_>((-PI/(float_)2) + eps, (PI/(float_)2) - eps, M); 
    //std::vector<float_> Angles = seq<float_>((-PI/(float_)2) + eps, (-PI/(float_)4) + eps, M); 
    //std::vector<float_> Heights = getHeights(M, eps, d);

    // Scattering map local definition
    //std::pair<std::vector<float_>, std::vector<float_>> P = getPointsInBranch_(KeyIndices[0], eps, d, M);
    //std::vector<float_> Heights = P.first;
    //std::vector<float_> Angles = P.second;
    //WriteVector<float_>("Data/BranchHeights.dat", Heights);
    //WriteVector<float_>("Data/BranchAngles.dat", Angles);
    
    // Scattering map local definition V2 (REJECTION METHOD)
    std::vector<int> TransmittingID, NTID, CustomID;
    std::vector<std::vector<int>> CustomIDs(4);
    if (d < 1) {
        TransmittingID = { 1,2,3,4,5,6,7,12,18 };
        NTID = { 0,8,9,10,11,13,14,15,16,17 };
        //CustomID = { 1,2,3,7,12,18 };
        //CustomID = { 1,2,12,18 };
	
	//CustomIDs[0] = {0,1,2,4,5,8,9,12,18};
	//CustomIDs[1] = {1,4,12};
	//CustomIDs[2] = {0,3,4,6,7,12,18};
	//CustomIDs[3] = {12};
	//CustomID = CustomIDs[KeyIndices[0]];

    } else {
        TransmittingID = { 4,5,6,12,18,19 };
        NTID = { 0,1,2,3,7,8,9,10,11,13,14,15,16,17 };

	//CustomIDs[0] = { 0,4,5,6,8,9,10,12,17,18,19 };
	//CustomIDs[1] = { 0,4,5,6,8,9,12,18,19 };
	//CustomIDs[2] = { 4,12,19 };
	//CustomIDs[3] = { 4,6,12,19 };
	//CustomID = CustomIDs[KeyIndices[0]];
    }
    //std::cout << "Running custom ID's with ID = " << CustomID << std::endl;

    //std::pair<std::vector<float_>, std::vector<float_>> P = getPointsInBranch_Rejection({KeyIndices[0]}, eps, d, M);
    //std::pair<std::vector<float_>, std::vector<float_>> P = getPointsInBranch_Rejection(TransmittingID, eps, d, M);
    ////std::pair<std::vector<float_>, std::vector<float_>> P = getPointsInBranch_Rejection(NTID, eps, d, M);
    //std::pair<std::vector<float_>, std::vector<float_>> P = getPointsInBranch_Rejection(CustomID, eps, d, M);
    //std::vector<float_> Heights = P.first;
    //std::vector<float_> Angles = P.second;
    //WriteVector<float_>("Data/BranchHeights.dat", Heights);
    //WriteVector<float_>("Data/BranchAngles.dat", Angles);
    // =============================== INITIAL CONDITIONS END ================================    
    

    // ================================================================================================
    // ============================= CONTINUOUS TIME SCATTERING MAP ==================================
    // ================================================================================================
    // Times to get. ~150 minutes for M = 1000 at 10^8 time units
    //float_ EndTime = 6e4; // This shouldn't be smaller than 1e5 (assuming total length of 1000)
    //float_ StartTime = EndTime / (float_)10;
    //std::vector<float_> Time = seq<float_>(StartTime, EndTime, 10);
    ////std::vector<float_> Time = { (float_)10000, (float_)25e4, (float_)1e5 };
    ////Time.insert(Time.begin(), (float_)1e4);
    //float_ dt = Time[3] - Time[2];
    ////std::vector<size_t> ttg = {200-1, 400-1, 600-1, 800-1, 1000-1};  // Direct index to @Time!!
    //std::vector<size_t> ttg = { 5, 9 };  // Direct index to @Time!!
    //std::cout << "t = " << Time[0] << ", " << Time[1] << ", ..., " << 
    //    Time.back() << ". dt = " << dt << std::endl;
    //std::cout << "ttg = " << index<float_>(Time, ttg) << std::endl;
    //
    //// Collecting moments
    //size_t nMoments = 8;
    //std::vector<Ensemble<float_>> Moments(nMoments, Ensemble<float_>(Time.size()));

    ////// Record initial conditions
    ////WriteVector<float_>(DataDir + "H0T0.dat", Heights);
    ////WriteVector<float_>(DataDir + "H0T0.dat", Angles);

    ////std::vector<float_> nu_FH_LB = { 0.247221, 0.276977, 0.327624, 0.410559, 0.420055, 0.425753, 0.446012, 0.492227, 0.498558, 0.523882, 0.541608, 0.563766, 0.6316};
    ////std::vector<float_> nu_FH_UB = { 0.250387, 0.280142, 0.330789, 0.413724, 0.423220, 0.428918, 0.449177, 0.495393, 0.501723, 0.527680, 0.546040, 0.567565, 0.6326}; 
    ////std::vector<float_> nu_IH_LB = { 0.330831, 0.348849, 0.359860, 0.389890, 0.436937, 0.483984, 0.500000, 0.630130, 0.664164, 0.685185, 0.694194, 0.725225 };  
    ////std::vector<float_> nu_IH_UB = { 0.335836, 0.353854, 0.364865, 0.394895, 0.441942, 0.488989, 0.505005, 0.635135, 0.669169, 0.690190, 0.699199, 0.730230 };

    //std::vector<float_> nu_FH_LB = { 0.498558 };
    //std::vector<float_> nu_FH_UB = { 0.501723 }; 
    //std::vector<float_> nu_IH_LB = { 0.664164 };  
    //std::vector<float_> nu_IH_UB = { 0.669169 };
    //Writer<float_> w_L(DataDir + "L.dat", true);
    //Writer<int_ll> w_R(DataDir + "R.dat", true);
    //Writer<float_> w_V(DataDir + "V.dat", true);
    //size_t skipCount = 0;

    //Random<float_> URHeight, URAngle;
    //// CONTINUOUS TIME SCATTERING MAP 
    //for (size_t k = 0; k < M; ++k) {
    //    if (10 * k % M == 0 && k != 0)
    //        std::cout << "Simulated " << k << " initial conditions." << std::endl;

    //    float_ Height = URHeight.getUniformRandom(eps, d - eps);
    //    float_ Angle = URAngle.getUniformRandom(-PI/(float_)2 + eps, (3*PI/(float_)2) - eps);

    //    //std::cout << Height << std::endl;
    //    ScatterTraj<float_> ST = getTrajectory(Height, Angle, d, Time);
    //    //std::cout << Height << std::endl;

    //    // It is likely skipping is relevant only if you want every time step
    //    if (ST.idx == Time.size()) {
    //        // Write data at every time step (only do small jobs)
    //        //WriteVector<float_>(DataDir + "H.dat", ST.Height);
    //        //WriteVector<int_ll>(DataDir + "R.dat", ST.Position);
    //        //WriteVector<float_>(DataDir + "A.dat", ST.Angle);
    //        //WriteVector<float_>(DataDir + "L.dat", ST.Tau);
    //        //WriteVector<int_ll>(DataDir + "BID.dat", ST.BranchID);
    //    
    //        // Write data at selected time steps
    //        std::vector<int_ll> R = index<int_ll>(ST.getPositions(), ttg);
    //        std::vector<float_> L = index<float_>(ST.getTaus(), ttg);

    //        //std::cout << R;
    //        //std::cout << L;

    //        std::vector<float_> Speed(ttg.size(), 0);
    //        for (size_t i = 0; i < ttg.size(); ++i) {
    //    	Speed[i] = (float_)R[i] / L[i];
    //        }

    //        //std::cout << Speed << std::endl;
    //        //if (all_gt<float_>(Speed, 1/(float_)5)) {
    //        if (abs(Speed[0]) > 1/(float_)5 && abs(Speed[1]) > 1/(float_)5) {
    //    	w_L.WriteRowVector(L);
    //    	w_R.WriteRowVector(R);
    //    	w_V.WriteRowVector(Speed);

    //    	//==========================================================================================
    //    	// Searching for trajectories which have a particular ballistic speed 
    //    	if (KeyIndices[2] == 0) { // Is the finite horizon configuration
    //    	    SearchAndWrite<float_>("Data/" + Key + "VB-FH_", Speed, nu_FH_LB, nu_FH_UB, 
    //    		    Height, Angle, d, false);
    //    	} else {
    //    	    SearchAndWrite<float_>("Data/" + Key + "VB-IH_", Speed, nu_IH_LB, nu_IH_UB, 
    //    		    Height, Angle, d, false);
    //    	}
    //    	//==========================================================================================
    //        }

    //        //std::vector<std::vector<float_>> DXp = getMomentPowers_INT<float_>(ST.getPositions(), nMoments);
    //        //for (size_t j = 0; j < Moments.size(); ++j) {
    //        //    Moments[j].Push(DXp[j]);
    //        //}
    //        skipCount++;
    //    }
    //}
    ////for (size_t i = 0; i < Moments.size(); ++i)
    ////    Moments[i].Write_Append(DataDir + "Moments.dat");
    //if (M - skipCount > 0)
    //    std::cout << M - skipCount << " trajectories were skipped." << std::endl;
    ////return 0;
    //// ================================================================================================
    //// ================================================================================================

    // ========================================================================
    // ===================== PERSISTENCE TIME SCATTERING MAP ==================
    // ========================================================================

    //Random<float_> URHeight, URAngle;
    //Logger<float_> Log(DataDir + "ST.dat", true);
    //for (size_t j = 0; j < Angles.size(); ++j) {
    //for (size_t k = 0; k < M; ++k) {
    //    if (10 * k % M == 0 && k != 0)
    //        std::cout << "Simulated " << k << " initial conditions." << std::endl;

    //    //float_ Height = URHeight.getUniformRandom(eps, d - eps);
    //    //float_ Angle = URAngle.getUniformRandom(-PI/(float_)2 + eps, (3*PI/(float_)2) - eps);

    //    float_ Height = Heights[k];
    //    float_ Angle = Angles[j];


    //    // Escape from a collection of regions
    //    //ScatterVar<float_> SV = getPersistenceState(Height, Angle, d, {KeyIndices[0]});
    //    //ScatterVar<float_> SV = getPersistenceState(Height, Angle, d, TransmittingID);
    //    //ScatterVar<float_> SV = getPersistenceState(Height, Angle, d, NTID);
    //    //ScatterVar<float_> SV = getPersistenceState(Height, Angle, d, CustomID);
    //    ScatterVar<float_> SV = JumpOnce(Height, Angle, d);
    //    //    Log.Stream << Height << " " << Angle << " " << SV.Height << " " << SV.Angle << " " << 
    //    //	SV.Position << " " << SV.Tau << " " << SV.BranchID << std::endl;
    //        Log.Stream << Height << " " << Angle * 180/PI << " " << SV.Tau << std::endl;
    //    //if (SV.Tau > 36000) {
    //    //    Log.Stream << SV.Height << " " << SV.Angle << " " << 
    //    //	SV.Position << " " << SV.Tau << " " << SV.BranchID << std::endl;
    //    //} else {
    //    //    M++;
    //    //}

    //    // Grid IC's2
    //    //int BID = whichBranch<float_>(Heights[k], Angles[j], d);
    //    //ScatterVar<float_> SV = getPersistenceState(Heights[k], Angles[j], d, { BID });

    //    //Log.Stream << Heights[k] << " " << Angles[j] << " " << 
    //    //	      SV.Position << " " << SV.Tau << " " << SV.BranchID << std::endl;
    //}
    ////Log.Stream << std::endl;
    //}

    // ========================================================================
    // ===================== DISCRETE TIME SCATTERING MAP =====================
    // ========================================================================

    // For plotting visuals - must use ScatterV2.h local definition
    Point<float_> IC = getPointInBranch_(KeyIndices[0], eps, d);
    ////Point<float_> IC = getPointInBranch_(18, eps, d);
    //Heights = { IC.x };
    //Angles = { IC.y };
    //Heights = { 1/(float_)5, 1/(float_)5 };
    //Angles = { -0.0712921, 0.0712921 };
    //Heights = { 0.15, 0.15 };
    //Angles = { -0.0555, 0.0555 };
    size_t N = 1;
    std::vector<size_t> itg = seq<size_t>(0, N, N+1);
    //std::vector<size_t> itg = {1};
    //std::vector<size_t> itg = seq<size_t>(0, N, N);

    //Logger<float_> Log_Persist(DataDir + "Persist.dat", true);
    //Logger<float_> Log(DataDir + "ST.dat", true);
    //Logger<float_> Log3d(DataDir + "H0T0RN.dat", true);
    Logger<float_> LogRegionCounts(DataDir + "RGC.dat", true);
    std::vector<int> TotalRegionsVisited(M, 0); // Unique regions visited along a trajectory with fixed angle and varying heights
    for (size_t j = 0; j < Angles.size(); ++j) {

    std::vector<int> Counts(20, 0);
    int TRV = 0;
    for (size_t k = 0; k < M; ++k) {
        //if (10 * k % M == 0 && k != 0)
        //    std::cout << "Simulated " << k << " initial conditions." << std::endl;

        //int BID = whichBranch<float_>(Heights[k], Angles[k], d);
        //ScatterVar<float_> SV = getPersistenceState(Heights[k], Angles[k], d, {BID});
        //Log_Persist.Stream << SV.Height << " " << SV.Angle << " " << 
        //	SV.Position << " " << SV.Tau << " " << SV.BranchID << std::endl;

        //ScatterTraj<float_> ST = getTrajectory(Heights[k], Angles[j], d, itg);
        ScatterTraj<float_> ST = getTrajectory(Heights[k], Angles[j], d, itg);
        //WriteColVector<int_ll>(DataDir + std::to_string(k) + "_R.dat", ST.getPositions());
        //WriteColVector<int>(DataDir + std::to_string(k) + "_BID.dat", ST.getBranchIDs());

        //Log.Stream << Heights[k] << " " << Angles[j] << " " << 
        //    ST.Position[1] << " " << ST.Tau[0] << " " << ST.BranchID[1] << std::endl;
        //for (size_t i = 0; i < itg.size(); ++i) {
        //    Log.Stream << ST.Height[i] << " " << ST.Angle[i] << " " << 
        //	ST.Position[i] << " " << ST.Tau[i] << " " << ST.BranchID[i] << std::endl;
        //}
        WriteVector<float_>(DataDir + "H.dat", ST.getHeights());
        WriteVector<int_ll>(DataDir + "R.dat", ST.getPositions());
        WriteVector<float_>(DataDir + "A.dat", ST.getAngles());
        WriteVector<float_>(DataDir + "L.dat", ST.getTaus());
        WriteVector<int>(DataDir + "BID.dat", ST.getBranchIDs());
        
        //std::vector<int_ll> pos = ST.getPositions();
        //std::vector<float_> ang = ST.getAngles();
        //std::vector<float_> hei = ST.getHeights();
        //Log3d.Stream << Heights[k] << " " << Angles[k] << " " <<  ang[1] << " " << hei[1] << " ";// << pos[1] << std::endl;
        //for (size_t i = 0; i < itg.size(); ++i)
        //    Log3d.Stream << pos[i] << " ";
        //Log3d.Stream << std::endl;

        std::vector<int> BIDs = ST.getBranchIDs(); // Regions visited by a trajectory with N collisions
        for (size_t i = 1; i < BIDs.size(); ++i) {
            Counts[BIDs[i]]++; // 
        }

        std::vector<int> uniqueBIDs = getUniqueElements(ST.getBranchIDs());
        //LogRegionCounts.Stream << Heights[k] << " " << Angles[j] << " " << uniqueBIDs.size()-1 << std::endl;

        TRV += (uniqueBIDs.size() - 1);
    }
        TotalRegionsVisited[j] = TRV; 
        //WriteVector<int>("Data/BID.dat", Counts);
    }
    WriteVector<int>("Data/TRV.dat", TotalRegionsVisited);
    // ========================= DISCRETE SCATTERING MAP END ===================================



    //// ========================================================
    //// Matching the scattering map with the billiard trajectory
    //// ========================================================
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";
    std::cout << RunTime << std::endl;
    
    //Simulation<float_> Sim {
    //    .Time = Time,
    //    .L = {0},
    //    .TimesToGet = {25000, 50000, 100000},
    //    .Symbols = {""},
    //    .DataDestination = "Data/" + Key,
    //    .RunTime = ".",
    //    .MomentStep = 2,
    //    .TotalTime = Time.back(),
    //    .nInitialConditions = M,
    //    .nCollisions = 0, // always reach max time before max coll
    //    .nMoments = 1};
    //getDisplacements_t<float_>(Table, Sim, ttg);

    //Particle<float_> P0 = Table.getInitialState(0, -1, 1e7);
    //float_ he = Heights[1];
    //float_ th = Angles[1];
    //P0.Pos = Point<float_>(0, he);
    //P0.Vel = Point<float_>(cos(th), sin(th));
    //Trajectory<float_> dTraj = Table.getTrajectory(P0, 2*N);
    //std::vector<float_> x = dTraj.getLiftedPositions_x();
    //std::vector<float_> y = dTraj.getLiftedPositions_y();
    ////std::vector<float_> fps = dTraj.getLiftedFreePaths();
    ////std::vector<float_> t(fps.size(), 0);
    ////std::partial_sum(fps.begin(), fps.end(), t.begin());
    ////WriteColVector<float_>(DataDir + "t.dat", t);
    //WriteColVector<float_>(DataDir + "x.dat", x);
    //WriteColVector<float_>(DataDir + "y.dat", y);
    ////std::vector<int> LocationFromTraj = getCrossingLocations(x, y);
    ////std::vector<float_> HeightsFromTraj = getCrossingHeights(x, y);
    ////WriteVector<float_>(DataDir + "RawH.dat", HeightsFromTraj);
    ////WriteVector<int>(DataDir + "RawX.dat", LocationFromTraj);
    // ====================================================

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
    

        // Correlation function (in the loop)
        //std::vector<float_> A = ST.getAngles();
        //for (size_t i = 0; i < CL.size(); ++i) {
        //    std::vector<float_> thetanm(N - CL.back() - CL[i], 0); // trailing zeros
        //    for (size_t j = 0; j < thetanm.size(); j++) 
        //        thetanm[j] = A[j] * A[j + CL[i]];
        //    THETAnm[i].Push(thetanm);
        //}

// Outside the loop
    //for (size_t i = 0; i < CL.size(); ++i) 
    //    THETAnm[i].Write(DataDir + "THETAn_" + std::to_string(i) + ".dat");

//
//		for (size_t i = 0; i < Vbals_FH.size(); ++i) {
//
//		    bool inRange = false;
//		    Logger<float_> w_VBLOGFH("Data/" + Key + "VB-FH_" + std::to_string(i) + ".dat", true);
//		    bool inRange = in_nh_ublb(abs(Speed[0]), Vbals_FH[i], Vbals_std_FH[i]) && 
//			in_nh_ublb(abs(Speed[1]), Vbals_FH[i], Vbals_std_FH[i]);
//
//		    Logger<float_> w_VBLOGIH("Data/" + Key + "VB-IH_" + std::to_string(i) + ".dat", true);
//		    bool inRange = in_nh_ublb(abs(Speed[0]), Vbals_IH[i], Vbals_std_IH[i]) && 
//		        in_nh_ublb(abs(Speed[1]), Vbals_IH[i], Vbals_std_IH[i]);
//
//		    if (inRange) { // Speed was inside range - now write trajectory data
//			//ScatterTraj<float_> STn = getTrajectory(Heights[k], Angles[k], d, itg);
//			// Write initial state and all other states
//			w_VBLOG.Stream << Heights[k] << " " << Angles[k] << " " << 
//			    0 << " " << 0 << " " << whichBranch<float_>(Heights[k], Angles[k], d) << std::endl;
//			for (size_t i = 1; i < itg.size(); ++i) {
//			    w_VBLOG.Stream << STn.Height[i] << " " << STn.Angle[i] << " " << 
//				STn.Position[i] << " " << STn.Tau[i] << " " << STn.BranchID[i] << std::endl;
//			}
//		    }
//		}
