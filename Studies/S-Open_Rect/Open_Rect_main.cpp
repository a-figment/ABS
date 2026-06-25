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
  // ======================= BOILER PLATE START ==========================
  std::cout << std::setprecision(std::numeric_limits<float_>::digits10) << std::showpoint;
  std::vector<int> KeyIndices;
  for (int i = 2; i < argc; ++i) 
      KeyIndices.push_back(std::atoi(argv[i]));
  DataKeys DK({"p"}, {10}, {"Partitions" });
  std::string Key = DK.MakeKey(KeyIndices);
  DK.Print();
  // ======================= BOILER PLATE END ==========================

  // ======================= BILLIARD DEFINITION START ==========================
  std::vector<std::unique_ptr<Obstacle<float_>>> Obstacles(1);


  Polygon<float_> Rect({0, 1, 1, 0}, {0, 0, 0.5, 0.5}, {1, 3}, false);

  Obstacles[0] = std::make_unique<Polygon<float_>>(Rect);
  Billiard<float_> Table(Obstacles);
  // ======================= BILLIARD DEFINITION END ==========================


  //std::vector<float_> Times = getExponentialTimes<float_>((float_)10000, 1000);
  std::vector<float_> Times = seq<float_>((float_)1, (float_)500000,  500000);

  // As many particles as possible to escape within the time limit
  Simulation<float_> Sim;
  Sim = Simulation<float_> {
      .Time = Times,
      .L = {0},
      .TimesToGet = { (float_)1000, (float_)100000 },
      .Symbols = {""},
      .DataDestination = "Data/" + Key,
      .RunTime = ".",
      .MomentStep = 1,
      .TotalTime = Times.back(),
      .nInitialConditions = 1000,
      .nCollisions = 0,
      .nMoments = 2
  };
  
 Sim.Print();
 //std::vector<size_t> Indices = getIndicesOfVectorValue({(float_)10, (float_)100, (float_)1000, (float_)10000}, Times);
 getDisplacementSummary_t(Table, Sim);
 //getPAC_t<float_>(Table, Sim);  
 //getDisplacements_t(Table, Sim, Indices);
 //WriteSampleTrajectory(Table, Sim);
 //ContinuousSystemSummary(Table, Sim);  
 //VelocityCorrelations<float_>(Table, Sim);  
 //CollisionSystemSummary(Table, Sim);  



  return 0;
}

