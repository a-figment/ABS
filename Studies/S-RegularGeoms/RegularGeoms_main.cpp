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
  DataKeys DK({"p", "d", "a"}, {1000, 2, 21}, {"Partitions", "[0] - FH, [1] - IF", "Alpha variation"});
  std::string Key = DK.MakeKey(KeyIndices);
  DK.Print();
  // ======================= BOILER PLATE END ==========================

  // ======================= BILLIARD DEFINITION START ==========================
  std::vector<float_> d = { -1 / ((float_)2), (float_)1 };
  std::vector<std::unique_ptr<Obstacle<float_>>> Obstacles(1);

  float_ PI = boost::math::constants::pi<float_>();

  // 21 values
  std::vector<float_> alpha = {PI / ((float_)12), PI / ((float_)11), PI / ((float_)10), PI / ((float_)9), PI / ((float_)8), PI / ((float_)7), PI / ((float_)6), PI / ((float_)5), PI / ((float_)4), PI / ((float_)3), PI / ((float_)2), (2*PI / ((float_)3)), (3*PI / ((float_)4)), (4*PI / ((float_)5)), (5*PI / ((float_)6)), (6*PI / ((float_)7)), (7*PI / ((float_)8)), (8*PI / ((float_)9)), (9*PI / ((float_)10)), (10*PI / ((float_)11)), (11*PI / ((float_)12)) };

  //Channel
  Polygon<float_> Channel = getParallelChannel<float_>(alpha[KeyIndices[2]], (1 / (float_)2) , d[KeyIndices[1]]);

  Obstacles[0] = std::make_unique<Polygon<float_>>(Channel);
  Billiard<float_> Table(Obstacles);
  // ======================= BILLIARD DEFINITION END ==========================

  Simulation<float_> Sim;
  Sim = Simulation<float_> {
      .Time = {0},
      .L = { 101 / ((float_)2) },
      .Symbols = {""},
      .DataDestination = "Data/" + Key,
      .RunTime = ".",
      .MomentStep = 0,
      .TotalTime = (float_)1e7,
      .nInitialConditions = 10000,
      .nCollisions = 0,
      .nMoments = 0};
  Sim.Write();
  FirstPassageTimes<float_>(Table, Sim);  



    //float_ T = (float_)1e6;
    //Simulation<float_> Sim {
    //    .Time = getExponentialTimes<float_>(T, 1000),
    //    .L = {0},
    //    .Symbols = {""},
    //    .DataDestination = "/fred/oz055/pathtohpc/ABS/S-RegularGeoms/" + Key,
    //    .RunTime = ".",
    //    .MomentStep = (float_)0.5,
    //    .TotalTime = T,
    //    .nInitialConditions = 10000,
    //    .nCollisions = 0,
    //    .nMoments = 20};
    //Sim.Write();
    // 
    //ContinuousSystemSummary(Table, Sim);  


  return 0;
}

