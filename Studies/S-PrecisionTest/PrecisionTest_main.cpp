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
  DataKeys DK({"p"}, {10}, {"Partitions"});
  std::string Key = DK.MakeKey(KeyIndices);
  DK.Print();
  // ======================= BOILER PLATE END ==========================

  // ======================= BILLIARD DEFINITION START ==========================
  std::vector<std::unique_ptr<Obstacle<float_>>> Obstacles(1);
  Polygon<float_> Channel = getParallelChannel<float_>(toRad<float_>(10), (1 / (float_)2) , (float_)0);

  Channel.Write("Data/BilliardTable.dat");
  Obstacles[0] = std::make_unique<Polygon<float_>>(Channel);
  Billiard<float_> Table(Obstacles);
  // ======================= BILLIARD DEFINITION END ==========================


  Simulation<float_> Sim;
  Sim = Simulation<float_> {
      .Time = getExponentialTimes<float_>((float_)1e6, 1000),
      .L = {0},
      .Symbols = {""},
      .DataDestination = "Data/" + Key,
      .RunTime = ".",
      .MomentStep = 1,
      .TotalTime = (float_)1e6,
      .nInitialConditions = 1000,
      .nCollisions = 0,
      .nMoments = 4};
  Sim.Write();
  
 ContinuousSystemSummary(Table, Sim);  

  return 0;
}

