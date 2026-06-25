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


std::vector<float_> RatioSelect(const float_& UB, const int& nRatios) {
    std::vector<float_> Ratio = seq<float_>(log2(1/(float_)4), log2(UB), nRatios);
    for (size_t i = 0; i < Ratio.size(); ++i)
	Ratio[i] = std::pow((float_)2, Ratio[i]); 
    return Ratio;
}


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
    DataKeys DK({"p", "a", "r"}, {1000, 100, 100}, {"PARTITIONS", "ALPHA_R = PI/4, PI/3, ..., 3PI/4", "RATIO"});
    std::string Key = DK.MakeKey(KeyIndices);
    DK.Print();
    // ======================= BOILER PLATE END ==========================
    


    // ======================= BILLIARD DEFINITION START ==========================
    std::vector<std::unique_ptr<Obstacle<float_>>> Obstacles(1);
    
    // Consts
    float_ PHI = boost::math::constants::phi<float_>();
    float_ PI = boost::math::constants::pi<float_>();

    // ================= ANGLES CONSIDERED ======================
    std::vector<float_> alpha_rational = { PI / (float_)4, PI / (float_)3, PI / (float_)2, 2*PI/(float_)3, 3*PI/(float_)4 };
    std::vector<float_> alpha(alpha_rational.size() * 2);
    for (size_t i = 0; i < alpha.size(); ++i) {
        if (i < (alpha_rational.size()))
            alpha[i] = alpha_rational[i];
        else
            alpha[i] = (alpha_rational[i - alpha_rational.size()] * (PHI - (3 / ((float_)(5)))));
    }
    float_ dy = (1/((float_)2)) * (cos(alpha[KeyIndices[1]] / (float_)2) / sin(alpha[KeyIndices[1]] / (float_)2)); 
    // ==========================================================

    // Different params require different ratio ranges
    std::vector<float_> ratio_FH = {  // Ratios for finite horizon configs
			1/(2*(1 - cos(alpha[KeyIndices[1]]))),
			1/(3*(1 - cos(alpha[KeyIndices[1]]))),
			1/(float_)5, 
			1/(float_)4, 
			7/(float_)24, 
			1/(float_)3, 
			9/(float_)20, 
			1/(float_)2, 
			2/(float_)3, 
			3/(float_)4, 
			(float_)1};


    // Special ratios for pi/3
    std::vector<float_> ratio_spec = { 19/(float_)10, 21/(float_)10, 39/(float_)10, 41/(float_)10, 59/(float_)10, 61/(float_)10, 79/(float_)10, 81/(float_)10 };
    std::vector<float_> ratio_IH = seq<float_>(3/(float_)2,10,18);
    std::vector<float_> ratio = ratio_IH;

    // Include all ratios into one vector and sort them
    ratio.insert(ratio.end(), ratio_FH.begin(), ratio_FH.end());
    ratio.insert(ratio.end(), ratio_spec.begin(), ratio_spec.end());
    std::sort(ratio.begin(), ratio.end());

    // Assign pore widths
    std::vector<float_> porewidths(ratio.size(), 0); 
    for (size_t i = 0; i < porewidths.size(); ++i) {
        porewidths[i] = ratio[i] * dy;
    }
    Polygon<float_> Channel = getParallelChannel2<float_>(dy, (1 / (float_)2), porewidths[KeyIndices[2]]);
    //Channel.Write("Data/" + Key + "BilliardTable.dat");
    Obstacles[0] = std::make_unique<Polygon<float_>>(Channel);
    Billiard<float_> Table(Obstacles);

    // ==================== SUMMARISE PARAMETER SELECTION ========================
    if ((size_t)KeyIndices[1] <= alpha_rational.size() - 1)  
        std::cout << "Running RATIONAL alpha = " + std::to_string(alpha[KeyIndices[1]] * 180/PI) << std::endl;
    else
        std::cout << "Running IRRATIONAL alpha = " + std::to_string(alpha[KeyIndices[1]] * 180/PI) << std::endl;

    std::cout << "Ratios = " << ratio << std::endl;
    std::cout << "Ratio selected = " << ratio[KeyIndices[2]] << std::endl;

    // ===========================================================================
    
    // ======================= BILLIARD DEFINITION END ==========================
    
    std::vector<float_> TIME = seq<float_>(0, 1e4, 11);
    //std::vector<float_> TIME = seq<float_>(0, 4e4, 11);
    std::vector<float_> ttg = TIME;
    ttg.erase(ttg.begin());

    Simulation<float_> Sim {
        .Time = TIME,
        .L = {0},
	//.TimesToGet = {(float_)16000, (float_)32000},
        .TimesToGet = {(float_)4000, (float_)8000},
        .Symbols = {""},
        .DataDestination = "Data/" + Key,
        .RunTime = ".",
        .MomentStep = 2,
        .TotalTime = TIME.back(),
        .nInitialConditions = 1000000,
        .nCollisions = 0, // always reach max time before max coll
        .nMoments = 1};
    

    getDisplacements_RS<float_>(Table, Sim, (float_)0);

    // Searching for trajectories belonging to ballistic fronts
    //std::vector<size_t> Indices = getIndicesOfVectorValue(Sim.TimesToGet, Sim.Time);
    //std::cout << "Collecting displacements at times t = " << index(Sim.Time, Indices) << " with indices = " << Indices << std::endl;

    ////float_ Vbal = 77/(float_)100;
    ////float_ Vbal = 75/(float_)100;
    //float_ Vbal = 50/(float_)100;
    ////float_ Vbal = 2/(float_)3;
    //float_ eps = 0.001;
    //for (size_t i = 0; i < Sim.nInitialConditions; ++i) {

    //    bool ThreshCond = true;
    //    Particle<float_> P0 = Table.getInitialState(0, -1, Sim.TotalTime);
    //    Particle<float_> Pn = P0;

    //    Trajectory<float_> Traj = Table.getTrajectory(P0, Sim.Time);

    //    // Rejection sample
    //    std::vector<float_> xDis = Traj.getDisplacements_x();
    //    float_ V = 0;
    //    for (size_t j = 0; j < Indices.size(); ++j) {
    //        float_ v = abs(xDis[Indices[j]]) / Sim.Time[Indices[j]]; // Make sure v = .5 at all times

    //        if (v < (Vbal - eps) || v > (Vbal + eps)) { 
    //    	ThreshCond = false;
    //    	break;
    //        } else {
    //    	V = v;
    //        }
    //    }

    //    // Write trajectory of the front
    //    //if (ThreshCond) {
    //    //    std::cout << i << "th particle had speed " << V << std::endl;
    //    //    Trajectory<float_> dTraj = Table.getTrajectory(Pn, 50000);
    //    //    WriteVector(Sim.DataDestination + "xPos.dat", dTraj.getLiftedPositions_x());
    //    //    WriteVector(Sim.DataDestination + "yPos.dat", dTraj.getLiftedPositions_y());
    //    //    std::vector<float_> fps = dTraj.getLiftedFreePaths();
    //    //    std::vector<float_> t(fps.size(), 0);
    //    //    std::partial_sum(fps.begin(), fps.end(), t.begin());
    //    //    WriteVector<float_>(Sim.DataDestination + "t.dat", t);
    //    //}
    //}





    //getDisplacementSummary_t<float_>(Table, Sim);

    //getMoments_t<float_>(Table, Sim);
    //WriteTrajectory_Perturb<float_>(Table, Sim);
    //getFront<float_>(Table, Sim);
    //getDisplacements_t<float_>(Table, Sim, { 2*1999, 2*3999, 2*7999, 2*15999 });
    //
    //float_ r = ratio[KeyIndices[2]];

    // PI3 FH
    //float_ vbal = 2 / (4 + 3*r);
    // PI2 IH
    //float_ vbal = 1 - (1 / (1 + r));


    //getFront<float_>(Table, Sim, 1/sqrt(3)); 

    //getfront<float_>(table, sim, vbal); // pi/3 for r14-r18 to r19
    //getFront<float_>(Table, Sim, vbal); // pi/2 for r11 to r19

    //getFront<float_>(Table, Sim, 0.555); // pi/2 for r11 to r19
    
    // PI2 FH
    //getFront<float_>(Table, Sim, 0.5020); // r1
    //getFront<float_>(Table, Sim, 0.50428); // r3
    //getFront<float_>(Table, Sim, 0.51451); // r5

    // PI3
    //getFront<float_>(Table, Sim, 0.735); // r16
    //getFront<float_>(Table, Sim, 0.915); // r18
    //getFront<float_>(Table, Sim, 0.401); // r19

    // PI4
    //getFront<float_>(Table, Sim, 0.2745); // r6
    //getFront<float_>(Table, Sim, 0.645); // r14
    //getFront<float_>(Table, Sim, 0.6625); // r18

    return 0;
}

