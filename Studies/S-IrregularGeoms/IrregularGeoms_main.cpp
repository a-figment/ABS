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
  DataKeys DK({"d", "a"}, {2, 500}, {"[0] - FH, [1] - IF", "Alpha variation"});
  std::string Key = DK.MakeKey(KeyIndices);
  DK.Print();
  // ======================= BOILER PLATE END ==========================

  // ======================= BILLIARD DEFINITION START ==========================
  std::vector<float_> d = { -1 / ((float_)2), (float_)1 };
  std::vector<std::unique_ptr<Obstacle<float_>>> Obstacles(1);

  float_ PI = boost::math::constants::pi<float_>();
  float_ PHI = boost::math::constants::phi<float_>();


  //
  std::vector<float_> alpha = seq<float_>(PI/((float_)12), 11*PI/((float_)12), 500);
  for (size_t i = 0; i < alpha.size(); ++i)
      alpha[i] = (PHI - 3/((float_)5)) * alpha[i];
  Polygon<float_> Channel = getParallelChannel<float_>(alpha[KeyIndices[1]], (1 / (float_)2) , d[KeyIndices[0]]);

  Obstacles[0] = std::make_unique<Polygon<float_>>(Channel);
  Billiard<float_> Table(Obstacles);
  // ======================= BILLIARD DEFINITION END ==========================

  Simulation<float_> Sim;
  Sim = Simulation<float_> {
      .Time = {0},
      .L = { 101/(float_)2 },
      .Symbols = {""},
      .DataDestination = "Data/" + Key,
      .RunTime = ".",
      .MomentStep = 0,
      .TotalTime = (float_)1e7,
      .nInitialConditions = 10000,
      .nCollisions = 0,
      .nMoments = 0};
  Sim.Write();
  FirstPassageTimes(Table, Sim);  


  return 0;
}

