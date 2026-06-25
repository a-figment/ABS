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
  DataKeys DK({"a", "d"}, {50, 50}, {"Internal angle (alpha)", "Width (d)"});
  std::string Key = DK.MakeKey(KeyIndices);
  DK.Print();
  // ======================= BOILER PLATE END ==========================

  // ======================= BILLIARD DEFINITION START ==========================
  std::vector<std::unique_ptr<Obstacle<float_>>> Obstacles(1);

  float_ PI = boost::math::constants::pi<float_>();
  float_ PHI = boost::math::constants::phi<float_>();


  // FIXED ALPHA, VARIABLE PORE WIDTH
  //std::vector<float_> alpha = { PI / ((float_)2), (PI / ((float_)2)) * (PHI - 3/((float_)5)) };
  //float_ dx = 1/((float_)2);
  //float_ dy = dx * (cos(alpha[KeyIndices[0]] / ((float_)2)) / (sin(alpha[KeyIndices[0]] / ((float_)2))));
  //std::vector<float_> widths = seq<float_>(1 / ((float_)10), dy, DK.DataKey.second[1]);

  //Polygon<float_> Channel = getParallelChannel2<float_>(dy, dx, widths[KeyIndices[1]]);

  //WriteVector<float_>("Data/Widths.dat", widths);

  // VARIABLE ALPHA, FIXED PORE WIDTH
  std::vector<float_> int_alpha; // 21 resonance values
  for (size_t i = 2; i <= 12; ++i)
      int_alpha.push_back(PI / ((float_)i));
  for (size_t i = 2; i <= 11; ++i)
      int_alpha.push_back(((float_)i * PI) / ((float_)i + (float_)1));
  std::vector<float_> alpha = seq<float_>(toRad<float_>(10), toRad<float_>(170), (size_t)DK.DataKey.second[1] - int_alpha.size());

  std::vector<float_> ALPHA = alpha;
  ALPHA.insert(ALPHA.end(), int_alpha.begin(), int_alpha.end());
  sort(ALPHA.begin(), ALPHA.end());

  Polygon<float_> Channel = getParallelChannel<float_>(ALPHA[KeyIndices[0]], (1 / (float_)2) , -1/((float_)2));
  WriteVector<float_>("Data/ALPHA.dat", ALPHA);


  Channel.Write("Data/" + Key + "BilliardTable.dat");
  Obstacles[0] = std::make_unique<Polygon<float_>>(Channel);
  Billiard<float_> Table(Obstacles);
  // ======================= BILLIARD DEFINITION END ==========================

  int M = 100000;
  Simulation<float_> Sim;
  Sim = Simulation<float_> {
      .Time = {0},
      .L = seq<float_>(3/((float_)2), (float_)50 + 1/((float_)2), 50),
      .Symbols = {""},
      .DataDestination = "Data/" + Key,
      .RunTime = ".",
      .MomentStep = 0,
      .TotalTime = (float_)1e3,
      .nInitialConditions = M,
      .nCollisions = 0,
      .nMoments = 0};
  Sim.Write();
  
  SFPT(Table, Sim);  

  return 0;
}

