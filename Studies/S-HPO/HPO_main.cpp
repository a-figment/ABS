#include <iostream>
#include <algorithm>
#include <chrono>
#include <thread>
#include "Billiard.h"
#include "Trajectory.h"
#include "OnlineStatistics.h"
#include "Misc.h"
#include "Helper.h"
#include "Ensemble.h"
#include "Simulation.h"
#include "DataKeys.h"
#include "Writer.h"
#include "Random.h"
#include "Logger.h"

// < single character prefix, max number of instantiations >
typedef float_ldub float_;
//typedef float_50 float_;
const int nDigits = std::numeric_limits<float_>::digits10;

//std::vector<float_> RatioSelect(const float_& UB, const int& nRatios) {
//    std::vector<float_> Ratio = seq<float_>(log2(1/(float_)4), log2(UB), nRatios);
//    for (size_t i = 0; i < Ratio.size(); ++i)
//	Ratio[i] = std::pow((float_)2, Ratio[i]); 
//    return Ratio;
//}

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
    DataKeys DK({"p", "a", "r", "t"}, {4000, 4000, 4000, 1000}, {"PARTITIONS", "ALPHA_R = PI/4, PI/3, ..., 3PI/4", "RATIO", "TIME"});
    std::string Key = DK.MakeKey(KeyIndices);
    DK.Print();
    // ======================= BOILER PLATE END ==========================
    
    // ======================= BILLIARD DEFINITION START ==========================
    std::vector<std::unique_ptr<Obstacle<float_>>> Obstacles(1);
    
    // Consts
    float_ PHI = boost::math::constants::phi<float_>();
    float_ PI = boost::math::constants::pi<float_>();

    // ================= ANGLES CONSIDERED ======================
    //std::vector<float_> alpha_rational = getRationalMultipleOfPi<float_>(7); 
    //std::vector<float_> alpha(alpha_rational.size() * 2);
    //for (size_t i = 0; i < alpha.size(); ++i) {
    //    if (i < (alpha_rational.size()))
    //        alpha[i] = alpha_rational[i];
    //    else
    //        alpha[i] = (alpha_rational[i - alpha_rational.size()] * (PHI - (3 / ((float_)(5)))));
    //}
    //std::vector<float_> alpha = seq<float_>(PI/(float_)30, 30*PI/(float_)31, 1000);
    //std::vector<float_> alpha = seq<float_>(PI/(float_)20, PI/(float_)2, 100);
    //std::vector<float_> alpha = seq<float_>(PI/(float_)2, 30*PI/(float_)31, 100);
    //std::vector<float_> alpha = seq<float_>(2*PI/(float_)3, 20*PI/(float_)21, 1000);
   

    std::vector<float_> alpha =  { PI/(float_)2 };
    
    float_ dx = (1/((float_)2));
    float_ dy = dx * cos(alpha[KeyIndices[1]] / (float_)2) / sin(alpha[KeyIndices[1]] / (float_)2);

    //std::vector<float_> ratio = seq<float_>(1/(float_)4, 10, 40);
    //std::vector<float_> ratio = seq<float_>(1/(float_)5, 3, 981);
    std::vector<float_> ratio = {1/(float_)4, 1/(float_)2, (float_)1, (float_)2};
    std::vector<float_> porewidths(ratio.size(), 0); //= { - 1/((float_)2), (float_)1 };
    for (size_t i = 0; i < porewidths.size(); ++i) {
        porewidths[i] = ratio[i] * dy;
    }
    float_ d = porewidths[KeyIndices[2]];

    //Polygon<float_> Channel = getParallelChannel2<float_>(0.3, 1, 0.1);
    Polygon<float_> Channel = getParallelChannel2<float_>(dy, (1 / (float_)2), d);
    //Polygon<float_> Channel = getParallelChannel<float_>(dy, (1 / (float_)2), porewidths[KeyIndices[2]]);
    //Channel.Write("Data/" + Key + "BilliardTable.dat");
    Obstacles[0] = std::make_unique<Polygon<float_>>(Channel);
    Billiard<float_> Table(Obstacles);

    // ==================== SUMMARISE PARAMETER SELECTION ========================
    std::cout << "alpha = " << alpha[KeyIndices[1]] * 180/PI << std::endl;
    //if (KeyIndices[1] <= alpha_rational.size() - 1)  
    //    std::cout << "alpha = " + alpha[KeyIndices[1]] * 180/PI << std::endl;
    //else
    //    std::cout << "Running IRRATIONAL alpha = " + alpha[KeyIndices[1]] * 180/PI << std::endl;

    std::cout << "Ratios = " << ratio;
    //for (size_t i = 0; i < 4; ++i)
    //    std::cout << ratio[i] << ", ";
    //std::cout << " ..., ";
    //for (size_t i = ratio.size()-3; i < ratio.size(); ++i)
    //    std::cout << ratio[i] << ", ";
    //std::cout << std::endl;
    // ===========================================================================
    
    // ======================= BILLIARD DEFINITION END ==========================
    std::vector<int> lengths = { 101, 1001, 10001, 100001 };
    std::vector<float_> TIME = seq<float_>(1, 1e4+1, lengths[KeyIndices[3]]);
    //std::vector<float_> ttg = TIME;
    //ttg.erase(ttg.begin());

    Simulation<float_> Sim {
        .Time = TIME,
        .L = {},
        .TimesToGet = {1},
        .Symbols = {""},
        .DataDestination = "Data/" + Key,
        .RunTime = ".",
        .MomentStep = 0,
        .TotalTime = TIME.back(),
        .nInitialConditions = 100000,
        .nCollisions = 100, // always reach max time before max coll
        .nMoments = 0};
    
    //getIHFlightTimes<float_>(Table, Sim);
    //getCollisionVelocity<float_>(Table, Sim);
    //SFPT<float_>(Table, Sim);
    //getHPO<float_>(Table, Sim); 
    //getNormalDistribution<float_>(10, 0, 3.14); 
    
    //WriteTrajectory_Perturb<float_>(Table, Sim);
    //WriteOne<float_>(Table, Sim, Point<float_>(0,0), Point<float_>(0,0), true);

    // Numerically getting the scattering angle
    float_ eps = 1e-8;
    Logger<float_> ScatterLog(Sim.DataDestination + "Scatter.dat", true);
    //Logger<float_> PosLog(Sim.DataDestination + "Pos.dat", true);

    Random<float_> H, A;
    std::vector<float_> H0 = seq<float_>(eps, d - eps, std::round(1000 * d));
    //H0 = { H0[100] };
    //std::vector<float_> H0 = {0.2};
    for (size_t i = 0; i < H0.size(); ++i) {
	Particle<float_> P0 = Table.getInitialState(0, -1, Sim.TotalTime);

	float_ Angle = A.getUniformRandom(-PI/(float_)2 + eps, (PI/(float_)2) - eps);
	float_ h0 = H.getUniformRandom(eps, d - eps);
	//float_ h0 = H0[i];
	P0.Pos = Point<float_>(-dx, h0);
	//P0.Vel = Point<float_>((float_)1, (float_)0); 
	P0.Vel = Point<float_>(cos(Angle), sin(Angle));
	P0.SegmentID = 5;

	std::vector<int> Itinerary;

	Point<float_> oldVel(0,0);
	size_t k = 0; // Number of collisions
	//PosLog.Stream << P0.Pos.x << " " << H0[i] << std::endl;
	while (abs(P0.LPos(Table.getDiameter()).x) <= 1/(float_)2) {
	    oldVel = P0.Vel;
	    Itinerary.push_back(P0.SegmentID);
	    Table.Update(P0);
	    k++;
	    //PosLog.Stream << P0.LPos(Table.getDiameter()).x << " " << P0.Pos.y << std::endl;
	}
	Itinerary.push_back(2);
	if (k == 0)
	    std::cout << "Did not update." << std::endl;
	//ScatterLog.Stream << h0 << " " << oldVel.x << " " << oldVel.y << " " << alpha[KeyIndices[1]] << std::endl;
	ScatterLog.Stream << Itinerary << std::endl;
    }




//    //std::vector<float_> h0 = { 0.25, 0.1 }; // pi/2, d = 0.3
//    //std::vector<float_> h0 = { 0.1, 0.8 }; // pi/3, d = 2.16506 (id = 230)
//    std::vector<float_> h0 = { 0.2825, 0.0556 }; // pi/3, d = 2.16506 (id = 230)
//    size_t M = h0.size();
//    Random<float_> H0;
//    for (size_t i = 0; i < M; ++i) {
//	Particle<float_> P0 = Table.getInitialState(0, -1, Sim.TotalTime);
//	P0.Pos = Point<float_>(-dx+eps, H0.getUniformRandom(eps, d - eps));
////	P0.Vel = Point<float_>((float_)1, (float_)0);
//	P0.Vel = Point<float_>((float_)0, (float_)1);
////	if (i >= 1)
////	    P0.Vel = Point<float_>((float_)0, (float_)1);
//
//	P0.SegmentID = 5;
//
//	P0.Pos.y = h0[i];
//
//	Trajectory<float_> dTraj = Table.getTrajectory(P0, 100000);
//	WriteVector<float_>(Sim.DataDestination + "xPos.dat", dTraj.getLiftedPositions_x());
//	WriteVector<float_>(Sim.DataDestination + "yPos.dat", dTraj.getLiftedPositions_y());
//    }


    //size_t M = 1;
    //std::vector<float_> Vbals(M, 0);
    //OnlineStatistics<float_> OS;
    //for (size_t i = 0; i < M; ++i) {
    //    Random<float_> H0;
    //    Particle<float_> P0 = Table.getInitialState(0, -1, Sim.TotalTime);
    //    P0.Pos = Point<float_>(-dx, H0.getUniformRandom(eps, d - eps));
    //    P0.Vel = Point<float_>((float_)1, (float_)0);
    //    P0.SegmentID = 5;

    //    P0.Pos.y = 0.3; // JHFSSL JHA:KJFDHA:WDIHWQPOIHIFEWU

    //    Trajectory<float_> dTraj = Table.getTrajectory(P, 100);
    //    Writevector<float_>(

    //    // Check first reflection vector
    //    Point<float_> iPos;
    //    Particle<float_> P_temp = P0;
    //    Table.UpdateCollisionState(P_temp);
    //    P_temp.Vel = ReflectSpecular<float_>(Channel.getInwardNormal(P_temp), P_temp.Vel);
    //    if (P_temp.Vel.x < 0) {
    //        //P0.Vel = Point<float_>(-1, 0);
    //        iPos = P_temp.Pos;
    //        //P0 = P_temp;
    //    }
    //        //P0.Vel = Point<float_>(cos(PI/(float_)3), sin(PI/(float_)3)); 

    //    Logger<float_> PosLog(Sim.DataDestination + "HPO.dat", true);
    //    PosLog.Stream << P0.Pos.x << " " << P0.Pos.y << " " << P0.SegmentID << " " << P0.TotalTime << std::endl;
    //    //while (P0.SegmentID != 2) {
    //    int n = 0;
    //    bool inPosyRange = false, inVelxRange = false, inVelyRange = false;

    //    do {
    //        //Table.UpdateCollisionState(P0);
    //        //P0.Vel = ReflectSpecular<float_>(Channel.getInwardNormal(P0), P0.Vel);
    //        
    //        Table.Update(P0);
    //        PosLog.Stream << P0.Pos.x << " " << P0.Pos.y << " " << P0.SegmentID << " " << P0.TotalTime << std::endl;
    //        n++;

    //        inPosyRange = in_nh(P0.Pos.y, iPos.y, eps);
    //        inVelxRange = in_nh(P0.Vel.x, -1/(float_)2, eps); 
    //        inVelyRange = in_nh(P0.Vel.y, sqrt(3)/(float_)2, eps); 
    //        std::cout << std::boolalpha << "n = " << n << ", " << inPosyRange << " " << inVelxRange << " " << inVelyRange << std::endl;
    //        std::cout << P0 << std::endl;

    //    } while (!(inPosyRange && inVelxRange && inVelyRange));

    //    Vbals[i] = 1 / P0.TotalTime;
    //    OS.Update(Vbals[i]);
    //}

    ////WriteVector(Sim.DataDestination + "Vbal.dat", Vbals);
    ////WriteVector(Sim.DataDestination + "Vbal.dat", OS.getMean());
    //Logger<float_> VBLog("Data/a" + std::to_string(KeyIndices[1]) + "-rX--" + "VBMean.dat", true);
    //VBLog.Stream << OS.getMean() << std::endl;




//============== COMPARISON OF JT ALGO WITH INFINITE HORIZON ==============
//    float_ lj_thresh = 500;
//
//    Writer<float_> _Xc(Sim.DataDestination + "Xc.dat", true);
//    Writer<float_> _Yc(Sim.DataDestination + "Yc.dat", true);
//    Writer<float_> _FP(Sim.DataDestination + "FP.dat", true);
//    Writer<float_> _Xt(Sim.DataDestination + "Xt.dat", true);
//    for (size_t i = 0; i < Sim.nInitialConditions; ++i) {
//    
//        Particle<float_> P0 = Table.getInitialState(0, -1, Sim.TotalTime);
//        Particle<float_> P1 = P0;
//    
//        Trajectory<float_> cTraj = Table.getTrajectory(P0, Sim.nCollisions);
//    
//	std::vector<float_> FP = cTraj.getLiftedFreePaths();
//	
//	bool has_large_jump = false;
//	for (size_t j = 0; j < FP.size(); ++j) {
//	    if (FP[j] > lj_thresh)  {
//		has_large_jump = true;
//		break;
//	    }
//	}
//
//	if (has_large_jump) {
//	    // Collisions
//            _Xc.WriteRowVector(cTraj.getLiftedPositions_x());
//            _Yc.WriteRowVector(cTraj.getLiftedPositions_y());
//            _FP.WriteRowVector(FP);
//    
//            // Continuous time
//            Trajectory<float_> conTraj = Table.getTrajectory(P1, Sim.Time);
//            _Xt.WriteRowVector(conTraj.getLiftedPositions_x());
//	}
//    
//    }

//============== TIME SERIES FOR DIFFERENT TIME STEPS  ==============
////Writer<float_> _Xt(Sim.DataDestination + "Xt.dat", true);
//Writer<float_> _Xc(Sim.DataDestination + "Xc.dat", true);
//Writer<float_> _Yc(Sim.DataDestination + "Yc.dat", true);
//Writer<float_> _FP(Sim.DataDestination + "FP.dat", true);
//for (size_t i = 0; i < Sim.nInitialConditions; ++i) {
//
//    Particle<float_> P0 = Table.getInitialState(0, -1, Sim.TotalTime);
//    Particle<float_> P1 = P0;
//
//    for (size_t j = 0; j < lengths.size(); ++j) {
//        Particle<float_> Psame = P0;
//        std::vector<float_> tim = seq<float_>(1, 1e7+1, lengths[j]);
//        Trajectory<float_> tTraj = Table.getTrajectory(Psame, tim);
//        WriteVector("Data/Positions/t" + std::to_string(j) + "-Xt.dat", tTraj.getLiftedPositions_x());
//    }
//
//    Trajectory<float_> cTraj = Table.getTrajectory(P0, Sim.nCollisions);
//    _Xc.WriteRowVector(cTraj.getLiftedPositions_x());
//    _Yc.WriteRowVector(cTraj.getLiftedPositions_y());
//    _FP.WriteRowVector(cTraj.getLiftedFreePaths());
//
//}
//===================================================================

//============== CHECKING MOMENTS WHEN INITIAL DISTRIBUTION IS NOT UNIFORM ==============
//    std::vector<float_> TIME = seq<float_>(1, 1e6+1, 1001);
//    int M = 100000;
//    int nm = 10;
//    std::string DataDir = "Data/ICVar/";
//
//    // Defining angles
//    Random<float_> Q1; // Quadrant 
//    std::vector<float_> init_vel(M, 0);
//    //float_ ba = PI / (float_)4;	    // primary angle
//    std::vector<float_> ba = { 0, PI / (float_)6, PI / (float_)4, PI / (float_)3 };
//    std::vector<float_> eps = { PI / (float_)180, 5 * PI / (float_)180 };   // epsilon either side of each angle
//
//
//
//    // CONSIDERING ALL COMBINATIONS OF 
//    for (size_t p = 0; p < ba.size(); ++p) {	  // Principle directions
//	for (size_t q = 0; q < eps.size(); ++q) { // Half-width of interval
//
//	    float_ shift = 0;
//	    for (int i = 0; i < M; ++i) {
//		if (i < M / 4) {
//		    shift = ba[p];
//		} else if (i >= M/4 && i < M/2) {
//		    shift = (ba[p] + PI/(float_)2);
//		} else if (i >= M/2 && i < 3*M/4) {
//		    shift = (ba[p] + PI);
//		} else if (i >= 3*M/4 && i < M) {
//		    shift = (ba[p] + 3*PI/(float_)2);
//		} else {
//		    std::cout << "Overflow" << std::endl;
//		}
//		init_vel[i] = Q1.getUniformRandom(shift - eps[q], shift + eps[q]);
//	    }
//	
//	    std::string prefix = "p" + std::to_string(p) + "-e" + std::to_string(q) + "-";
//	    WriteVector<float_>(DataDir + prefix + "iv.dat", init_vel);
//	
//	
//	    // For the moments
//	    std::vector<Ensemble<float_>> Moments(nm, Ensemble<float_>(TIME.size()));
//	    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
//	    for (int i = 0; i < M; ++i)
//	    {
//	        if (10 * i % M == 0 && i != 0)
//	            std::cout << "Simulated " << i << " initial conditions." << std::endl;
//	        Particle<float_> P0 = Table.getInitialState(0, -1, TIME.back());
//		P0.Vel = Point<float_>(cos(init_vel[i]), sin(init_vel[i]));
//	        Trajectory<float_> Traj = Table.getTrajectory(P0, TIME);
//	        std::vector<std::vector<float_>> DXp = Traj.getDisplacements_x(1/(float_)2, nm);
//	        for (int j = 0; j < nm; ++j)
//	            Moments[j].Push(DXp[j]);
//	    }
//	    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
//	    std::cout << "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n" << std::endl;
//	
//	    for (size_t i = 0; i < Moments.size(); ++i)
//	        Moments[i].Write(DataDir + prefix + "Moment_" + std::to_string(i) + ".dat");
//	
//	}
//    }



    return 0;
}

