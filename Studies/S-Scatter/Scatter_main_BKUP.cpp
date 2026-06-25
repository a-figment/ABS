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


// @x,y 1x2 vectors.
template<typename T>
T getCrossingHeight(const std::vector<T>& x, const std::vector<T>& y) {
    std::vector<int> ix = intf(x);
    int mid = 0;
    T yint = 0;
    if (ix[0] != ix[1]) {

	T m = (y[1] - y[0])/(x[1] - x[0]);

	if ((ix[0] == -1 && ix[1] == 1) || (ix[0] == 1 && ix[1] == -1)) {
	    mid = 0;
	} else if (ix[0] > 0 && ix[1] > 0) {
	    mid = std::min(ix[0], ix[1]);
	} else if (ix[0] < 0 && ix[1] < 0) {
	    mid = std::max(ix[0], ix[1]);
	}
	yint = m * ((T)mid - x[0]) + y[0];
	return yint;
    } else {
	return INFINITY;
    }
}

// @x,y 1x2 vectors.
template<typename T>
int getCrossingLocation(const std::vector<T>& x, const std::vector<T>& y) {
    std::vector<int> ix = intf(x);
    int mid = 0;
    T yint = 0;
    if (ix[0] != ix[1]) {

	T m = (y[1] - y[0])/(x[1] - x[0]);

	if ((ix[0] == -1 && ix[1] == 1) || (ix[0] == 1 && ix[1] == -1)) {
	    mid = 0;
	} else if (ix[0] > 0 && ix[1] > 0) {
	    mid = std::min(ix[0], ix[1]);
	} else if (ix[0] < 0 && ix[1] < 0) {
	    mid = std::max(ix[0], ix[1]);
	}
	yint = m * ((T)mid - x[0]) + y[0];
	return mid;
    } else {
	return std::numeric_limits<int>::max();
    }
}

template<typename T>
std::vector<T> getCrossingHeights(const std::vector<T>& x, const std::vector<T>& y) {

    std::vector<T> Heights;
    Heights.reserve(x.size());
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
    for (size_t i = 0; i < x.size()-1; ++i) {
	std::vector<T> X = {x[i], x[i+1]};
	std::vector<T> Y = {y[i], y[i+1]};
	int F = getCrossingLocation(X, Y);
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
void Recursion(const T& d, const T& eps, const size_t& M,  const std::vector<T>& ttg) {
    Random<T> UR;
    for (size_t k = 0; k < M; ++k) {

	float_ Angle = UR.getUniformRandom(0, 2*PI - eps);
	float_ Height = UR.getUniformRandom((float_)0 + eps, d - eps);

	// Initial conditions
	float_ h0 = Height;
	float_ T0 = Angle;
	float_ tau = 0;
	int_ll x0 = 0;
	int bid = -1;
	bool sym = false;

	size_t tidx = 0;
	float_ TAU = tau;

	Scatterer<float_> S(d, h0, T0, x0, bid);
	while (TAU < ttg.back()) {
	    S.Evolve(h0, T0, x0, tau, bid, sym);
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
	Angles[i] = UR.getUniformRandom(0, 2*PI - eps);
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
    T T0 = URA.getUniformRandom(0, 2*PI - eps);
    T tau = 0;
    int_ll x0 = 0;
    int bid = -1;
    bool sym = false;

    T foundH0 = -1;
    T foundT0 = 0;
    do {
	h0 = URH.getUniformRandom(eps, d - eps);
	//T0 = URA.getUniformRandom((3*PI/(T)2) + eps, 2*PI);
	//T0 = URA.getUniformRandom((PI/(T)2) + eps, (3*PI/(T)2) - eps);
	T0 = URA.getUniformRandom(0, 2*PI - eps);
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
void JumpTimes(const size_t& M, const T& minTime, const std::string& DataDir) {
    
    T d = 1/(T)2;
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
    float_ d = 1/(float_)2;
    Polygon<float_> Channel({0, dx, 2*dx, 2*dx, dx, 0}, {0, dy, 0, d, d + dy, d}, {2, 5}, false);
    std::cout << Channel << std::endl;

    Obstacles[0] = std::make_unique<Polygon<float_>>(Channel);
    Billiard<float_> Table(Obstacles);

    // =============================== SIMULATION START ===============================    
    std::string DataDir = "Data/" + Key;
    //std::string DataDir = "/fred/oz055/jorchard_data/ABS/Scatter/" + Key;

    float_ eps = 1e-12;
    size_t M = 1000000;

    
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

    // Uniform random IC's (the default)
    std::vector<float_> Angles = getAngles(M, eps);
    std::vector<float_> Heights = getHeights(M, eps, d);

    // User defined IC
    //std::vector<float_> Heights = { 7/(float_)100 };
    //std::vector<float_> Angles(M, (16*PI/(float_)9)); // = { atan(1/(float_)2) };
    //std::vector<float_> Heights = getHeights(M, eps, d);

    // Scattering map local definition
    //std::pair<std::vector<float_>, std::vector<float_>> P = getPointsInBranch(KeyIndices[0], eps, d, M);
    //std::vector<float_> Heights = P.first;
    //std::vector<float_> Angles = P.second;

    //std::pair<std::vector<float_>, std::vector<float_>> P = getPointsInBranch(KeyIndices[0], eps, d, 10000, M);
    //std::vector<float_> Heights = P.first;
    //std::vector<float_> Angles = P.second;
    // =============================== INITIAL CONDITIONS END ================================    
    
    //return 0;
    //JumpTimes<float_>(M, (float_)0, DataDir);
    //return 0;


    // ==========================================================================
    // ===================== CONTINUOUS TIME SCATTERING MAP =====================
    // ==========================================================================
    // Times to get. ~150 minutes for M = 1000 at 10^8 time units
    float_ EndTime = 1e8; // This shouldn't be smaller than 1e5 (assuming total length of 1000)
    float_ StartTime = EndTime / (float_)1000;
    std::vector<float_> Time = seq<float_>(StartTime, EndTime, 1000);
    Time.insert(Time.begin(), (float_)1e4);
    float_ dt = Time[1] - Time[0];
    std::vector<size_t> ttg = {0, 10, 100, 200, 400, 800, 1000}; 
    
    // Collecting moments
    size_t nMoments = 8;
    std::vector<Ensemble<float_>> Moments(nMoments, Ensemble<float_>(Time.size()));

    // Record initial conditions
    //WriteVector<float_>(DataDir + "H0T0.dat", Heights);
    //WriteVector<float_>(DataDir + "H0T0.dat", Angles);

    size_t skipCount = 0;

    // CONTINUOUS TIME SCATTERING MAP 
    //for (size_t k = 0; k < M; ++k) {
    //    if (10 * k % M == 0 && k != 0)
    //        std::cout << "Simulated " << k << " initial conditions." << std::endl;

    //    // IC's
    //    float_ h0 = Heights[k];
    //    float_ T0 = Angles[k];
    //    float_ tau = 0;
    //    int_ll x0 = 0;
    //    int bid = -1;
    //    bool sym = false;

    //    // Evolve in continuous time
    //    size_t tidx = 0;
    //    float_ TAU = tau;  // Current time
    //    float_ pTAU = TAU; // Previous time

    //    Scatterer<float_> S(d, h0, T0, x0, bid);	// Evolves the map
    //    ScatterVar<float_> SV(h0, T0, pTAU, x0, bid);      // Variables of the map
    //    ScatterTraj<float_> ST(Time.size()); // Stores variables 
    //    bool skip = false;

    //    // Evolve map until time exceeds target time
    //    while (TAU < Time.back()) {

    //        S.Evolve_(h0, T0, x0, tau, bid, sym);
    //        sym = false;
    //        TAU += tau;
    //        
    //        // 
    //        if (tau > 2*dt) {
    //    	std::cout << "Crossing time tau = " << tau << " was longer than step time = " << dt << std::endl;
    //    	std::cout << S << std::endl;
    //    	SV.Print();
    //    	skip = true;
    //    	break;
    //        }

    //        if (tidx < Time.size() && TAU >= Time[tidx]) { // Exceeded a target time
    //    	// Is current iterate closer to Time than previous iterate?
    //    	if (abs(TAU - Time[tidx]) < abs(pTAU - Time[tidx])) {
    //    	    ST.Update(h0, T0, TAU, x0, bid);
    //    	} else {
    //    	    ST.Update(SV);
    //    	}
    //            tidx++;
    //        }
    //        pTAU = TAU;
    //        SV.Set(h0, T0, pTAU, x0, bid);
    //    }

    //    if (!skip) {
    //        //WriteVector<float_>(DataDir + "H.dat", ST.Height);
    //        //WriteVector<int_ll>(DataDir + "R.dat", ST.Position);
    //        //WriteVector<float_>(DataDir + "A.dat", ST.Angle);
    //        //WriteVector<float_>(DataDir + "L.dat", ST.Tau);
    //        //WriteVector<int_ll>(DataDir + "BID.dat", ST.BranchID);
    //        //std::vector<int_ll> pos = ST.Position;
    //    
    //        std::vector<std::vector<float_>> DXp = getMomentPowers_INT<float_>(ST.getPositions(), nMoments);
    //        for (size_t j = 0; j < Moments.size(); ++j) {
    //            Moments[j].Push(DXp[j]);
    //        }
    //        WriteVector<int_ll>(DataDir + "R.dat", index<int_ll>(ST.getPositions(), ttg));
    //        WriteVector<float_>(DataDir + "L.dat", index<float_>(ST.getTaus(), ttg));
    //        WriteVector<float_>(DataDir + "A.dat", index<float_>(ST.getAngles(), ttg));
    //        skipCount++;
    //    }
    //}
    //for (size_t i = 0; i < Moments.size(); ++i)
    //    Moments[i].Write_Append(DataDir + "Moments.dat");
    //std::cout << M - skipCount << " trajectories were skipped." << std::endl;



    // ========================================================================
    // ===================== DISCRETE TIME SCATTERING MAP =====================
    // ========================================================================

    //Point<float_> IC = getPointInBranch(KeyIndices[0], eps, d);
    //Point<float_> IC = Point<float_>(0.236067369601369095, 32.5821945296272611*PI/(float_)180);
    //std::cout << IC.x << " " << IC.y * 180/PI << std::endl;

    // Correlations in the departure angles
    size_t N = 10;
    std::vector<size_t> itg = seq<size_t>(0, N, N+1);
    //std::vector<size_t> itg = seq<size_t>(1, N, N);
    
    // Angle-angle correlations
    //std::vector<size_t> CL = { 1, 2, 3, 4, 5, 10, 20, 25 }; // Correlation lengths 
    //std::vector<size_t> CL = seq<size_t>(1, 100, 100);
    //assert(CL.back() < N);
    //std::vector<Ensemble<float_>> THETAnm(CL.size(), Ensemble<float_>(N - CL.back()));

    //Logger<float_> Log(DataDir + "BHA.dat", true);
    for (size_t k = 0; k < M; ++k) {
        if (10 * k % M == 0 && k != 0)
            std::cout << "Simulated " << k << " initial conditions." << std::endl;

        // IC's
        float_ h0 = Heights[k];
        float_ T0 = Angles[k];
        float_ tau = 0;
        int_ll x0 = 0;
        int bid = -1;
        bool sym = false;

        size_t j = 0;
        float_ TAU = 0;

	Scatterer<float_> S(d, h0, T0, x0, bid);	// Evolves the map

	ScatterTraj<float_> ST(itg.size()-1);
	for (size_t i = 0; i < itg.back(); ++i) {
	    if (j < itg.size() && itg[j] == i) {
		ST.Update(h0, T0, TAU, x0, bid); // 
		j++;
	    }
	    S.Evolve_(h0, T0, x0, tau, bid, sym);
	    TAU += tau;
	    sym = false; // MUST RESET!!!!
	}


	// Correlation function
	//std::vector<float_> A = ST.getAngles();
	//for (size_t i = 0; i < CL.size(); ++i) {
        //    std::vector<float_> thetanm(N - CL.back() - CL[i], 0); // trailing zeros
	//    for (size_t j = 0; j < thetanm.size(); j++) 
        //        thetanm[j] = A[j] * A[j + CL[i]];
        //    THETAnm[i].Push(thetanm);
	//}


        //WriteColVector<int_ll>(DataDir + std::to_string(k) + "_R.dat", ST.getPositions());
        //WriteColVector<int>(DataDir + std::to_string(k) + "_BID.dat", ST.getBranchIDs());

        //WriteVector<float_>(DataDir + "H.dat", ST.getHeights());
        //WriteVector<int_ll>(DataDir + "R.dat", ST.getPositions());
        //WriteVector<float_>(DataDir + "A.dat", ST.getAngles());
        //WriteVector<float_>(DataDir + "L.dat", ST.getTaus());
        //WriteVector<int>(DataDir + "BID.dat", ST.getBranchIDs());

	//std::vector<int> uniqueBIDs = getUniqueElements(ST.getBranchIDs());
	//if (uniqueBIDs.size() - 1 < 7)
	    //Log.Stream << Heights[k] << " " << Angles[k] * 180/PI << " " << uniqueBIDs.size()-1 << std::endl;
    }
    for (size_t i = 0; i < CL.size(); ++i) 
        THETAnm[i].Write(DataDir + "THETAn_" + std::to_string(i) + ".dat");

    // ========================= DISCRETE SCATTERING MAP END ===================================



    // ========================================================
    // Matching the scattering map with the billiard trajectory
    // ========================================================
    Particle<float_> P0 = Table.getInitialState(0, -1, 1e7);
    float_ he = Heights[0];
    float_ th = Angles[0];
    P0.Pos = Point<float_>(0, he);
    P0.Vel = Point<float_>(cos(th), sin(th));

    Trajectory<float_> dTraj = Table.getTrajectory(P0, 1000);
    std::vector<float_> x = dTraj.getLiftedPositions_x();
    std::vector<float_> y = dTraj.getLiftedPositions_y();
    std::vector<float_> fps = dTraj.getLiftedFreePaths();
    std::vector<float_> t(fps.size(), 0);
    std::partial_sum(fps.begin(), fps.end(), t.begin());
    //WriteColVector<float_>(DataDir + "t.dat", t);
    WriteColVector<float_>(DataDir + "x.dat", x);
    WriteColVector<float_>(DataDir + "y.dat", y);

    std::vector<int> LocationFromTraj = getCrossingLocations(x, y);
    std::vector<float_> HeightsFromTraj = getCrossingHeights(x, y);
    WriteVector<float_>(DataDir + "RawH.dat", HeightsFromTraj);
    WriteVector<int>(DataDir + "RawX.dat", LocationFromTraj);
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

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";
    std::cout << RunTime << std::endl;
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
    

