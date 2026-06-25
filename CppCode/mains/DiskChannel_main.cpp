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
    DataKeys DK({"p", "a", "r", "s"}, {1000, 10, 100, 10}, {"PARTITIONS", "ALPHA_R = PI/4, PI/3, ..., 3PI/4", "RATIO", "RADIUS"});
    std::string Key = DK.MakeKey(KeyIndices);
    DK.Print();
    // ======================= BOILER PLATE END ==========================
    


    // ======================= BILLIARD DEFINITION START ==========================
    
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
    int nRatios = 20;
    std::vector<float_> ratio(nRatios);
    size_t a = KeyIndices[1];
    if (a == 0 || a == (0 + alpha_rational.size())) {
	  ratio = RatioSelect(10, nRatios);
    } else if (a == 1 || a == (1 + alpha_rational.size())) {
	  ratio = RatioSelect(2, nRatios);
    } else if (a == 2 || a == (2 + alpha_rational.size())) {
	  ratio = RatioSelect(10, nRatios);
    } else if (a == 3 || a == (3 + alpha_rational.size())) {
	  ratio = RatioSelect(10, nRatios);
    } else if (a == 4 || a == (4 + alpha_rational.size())) {
	  ratio = RatioSelect(10, nRatios);
    } else if (a == 5 || a == (5 + alpha_rational.size())) {
	  ratio = RatioSelect(10, nRatios);
    } else {
	  std::cout << "Invalid option. No ratios defined!" << std::endl;
    }

    std::vector<float_> porewidths(ratio.size(), 0); //= { - 1/((float_)2), (float_)1 };
    for (int i = 0; i < porewidths.size(); ++i) {
        porewidths[i] = ratio[i] * dy;
    }

    std::vector<float_> Radii = { (float_)0, 1/(float_)10 };
    float_ DiskRadius = Radii[KeyIndices[3]];

    float_ ca2 = cos(alpha[KeyIndices[2]] / (float_)2);
    float_ sa2 = sin(alpha[KeyIndices[2]] / (float_)2);

    Polygon<float_> TC = getParallelChannel2<float_>(dy, (1 / (float_)2), porewidths[KeyIndices[2]]); // Temp channel
    float_ dx = TC.getVertices()[2].x;
    float_ d = TC.getVertices()[3].y;


    // -------- RECONSTRUCT CHANNEL --------
    float_ H  = TC.getVertices()[3].y + (DiskRadius / sa2);  // Left and right distance to center
    float_ Hc = TC.getVertices()[1].y - (DiskRadius / sa2);  // Center distance to center

    // Center circles points of tangency
    Point<float_> CLp(-DiskRadius * ca2, DiskRadius * sa2 + Hc);   // Center, left tangent point
    Point<float_> CRp( DiskRadius * ca2, DiskRadius * sa2 + Hc);   // Center, right tangent point

    // Left circle (RIGHT) point of tangency
    Point<float_> LRp( DiskRadius * ca2 - TC.getVertices()[0].x, - DiskRadius * sa2 + H);

    // Right circle (LEFT) point of tangency
    Point<float_> RLp(-DiskRadius * ca2 + TC.getVertices()[0].x, - DiskRadius * ca2 + H);

    float_ theta = atan(TC.getVertices()[1].y / TC.getVertices()[2].x);

    std::vector<float_> xv = {-dx,       CLp.x, CRp.x, dx,        dx,             RLp.x, TC.getVertices()[1].x,            RLp.x, -dx,             -dx      };
    std::vector<float_> yv = {(float_)0, CLp.y, CRp.y, (float_)0, H - DiskRadius, RLp.y, dy + d - DiskRadius * cos(theta), RLp.y, H - DiskRadius, (float_)0};

    std::cout << xv.size() << " -- " << yv.size() << std::endl;
    Polygon<float_> Channel(xv, yv, {3,6}, false, true);
    Channel.Write("Data/ModifiedChannel.dat");

    // Left, center and right circles
    Point<float_> Lc(TC.getVertices()[0].x, H);
    Point<float_> Cc(TC.getVertices()[1].x, Hc);
    Point<float_> Rc(TC.getVertices()[2].x, H);
    Circle<float_> LeftCircle(DiskRadius,Lc);
    Circle<float_> CenterCircle(DiskRadius,Cc);
    Circle<float_> RightCircle(DiskRadius,Rc);


    std::vector<std::unique_ptr<Obstacle<float_>>> Obstacles(4); // Channel, center and left/right circles
    Obstacles[0] = std::make_unique<Polygon<float_>>(Channel);
    Obstacles[1] = std::make_unique<Circle<float_>>(CenterCircle);
    Obstacles[2] = std::make_unique<Circle<float_>>(LeftCircle);
    Obstacles[3] = std::make_unique<Circle<float_>>(RightCircle);
    Billiard<float_> Table(Obstacles);

    // ==================== SUMMARISE PARAMETER SELECTION ========================
    if (KeyIndices[1] <= alpha_rational.size() - 1)  
        std::cout << "Running RATIONAL alpha = " + std::to_string(alpha[KeyIndices[1]] * 180/PI) << std::endl;
    else
        std::cout << "Running IRRATIONAL alpha = " + std::to_string(alpha[KeyIndices[1]] * 180/PI) << std::endl;

    std::cout << "Ratios = ";
    for (size_t i = 0; i < 4; ++i)
	std::cout << ratio[i] << ", ";
    std::cout << " ..., ";
    for (size_t i = ratio.size()-3; i < ratio.size(); ++i)
	std::cout << ratio[i] << ", ";
    std::cout << std::endl;

    std::cout << "Simulating disk of radius s = " << DiskRadius << "." << std::endl;
    // ===========================================================================
    
    // ======================= BILLIARD DEFINITION END ==========================
    
    std::vector<float_> TIME = { (float_)0, (float_)4000, (float_)8000, (float_)16000, (float_)32000, (float_)64000 };


    Simulation<float_> Sim {
        .Time = TIME,
        .L = {0},
        .TimesToGet = {0},
        .Symbols = {""},
        //.DataDestination = "/fred/oz055/jorchard_data/ABS/S-BallisticFront/" + Key,
        .DataDestination = "Data/" + Key,
        .RunTime = ".",
        .MomentStep = 0,
        .TotalTime = TIME.back(),
        .nInitialConditions = 10000,
        .nCollisions = 0, // always reach max time before max coll
        .nMoments = 0};
    //Sim.Write();
    
    getDisplacements_t<float_>(Table, Sim, { 2*1999, 2*3999, 2*7999, 2*15999 });//,6});  

   return 0;
}

