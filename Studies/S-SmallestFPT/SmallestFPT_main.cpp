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
  DataKeys DK({"a"}, {521}, {"Alpha variation"});
  std::string Key = DK.MakeKey(KeyIndices);
  DK.Print();
  // ======================= BOILER PLATE END ==========================

  // ======================= BILLIARD DEFINITION START ==========================
  std::vector<std::unique_ptr<Obstacle<float_>>> Obstacles(1);

  float_ PI = boost::math::constants::pi<float_>();
  float_ PHI = boost::math::constants::phi<float_>();

  //std::vector<float_> alpha = seq<float_>(PI / ((float_)6), 5 * PI / ((float_)6), DK.DataKey.second[1]);
  std::vector<float_> alpha = seq<float_>(toRad<float_>(10), toRad<float_>(170), 500);
  //std::vector<float_> alpha = pseq<float_>(toRad<float_>(PHI * 10), toRad<float_>(100 * PHI), DK.DataKey.second[1]);
 

  // 21 values
  std::vector<float_> resonance_values = {PI / ((float_)12), PI / ((float_)11), PI / ((float_)10), PI / ((float_)9), PI / ((float_)8), PI / ((float_)7), PI / ((float_)6), PI / ((float_)5), PI / ((float_)4), PI / ((float_)3), PI / ((float_)2), (2*PI / ((float_)3)), (3*PI / ((float_)4)), (4*PI / ((float_)5)), (5*PI / ((float_)6)), (6*PI / ((float_)7)), (7*PI / ((float_)8)), (8*PI / ((float_)9)), (9*PI / ((float_)10)), (10*PI / ((float_)11)), (11*PI / ((float_)12)) };

  std::vector<float_> ALPHA = alpha;
  ALPHA.insert(ALPHA.end(), resonance_values.begin(), resonance_values.end());
  sort(ALPHA.begin(), ALPHA.end());

  WriteVector<float_>("Data/ALPHA.dat", ALPHA);

  // FINITE HORIZON
  Polygon<float_> Channel = getParallelChannel<float_>(ALPHA[KeyIndices[0]], (1 / (float_)2) , -1/((float_)2));

  // INFINITE HORIZON
  //Polygon<float_> Channel = getParallelChannel<float_>(ALPHA[KeyIndices[0]], (1 / (float_)2) , (float_)1);

  //Channel.Write("Data/" + Key + "BilliardTable.dat");
  Obstacles[0] = std::make_unique<Polygon<float_>>(Channel);
  Billiard<float_> Table(Obstacles);
  // ======================= BILLIARD DEFINITION END ==========================

  int M = 1000000;

  float_ H = 1/((float_)2);

  Simulation<float_> Sim;
  Sim = Simulation<float_> {
      .Time = {0},
      .L = { (float_)50 + H, (float_)100 + H, (float_)250 + H },
      .Symbols = {""},
      .DataDestination = "Data/" + Key,
      .RunTime = ".",
      .MomentStep = 0,
      .TotalTime = (float_)5000,
      .nInitialConditions = M,
      .nCollisions = 0,
      .nMoments = 0};
  Sim.Write();
  
  SFPT<float_>(Table, Sim);  

  return 0;
}

