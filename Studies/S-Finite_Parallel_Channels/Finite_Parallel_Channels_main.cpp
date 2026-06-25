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


// < single character prefix, max number of instantiations >
typedef float_ldub float_;
const int nDigits = std::numeric_limits<float_>::digits10;


struct DataKeys {
    
    DataKeys(const std::vector<std::string>& s, const std::vector<uint>& maxs) {
        assert(s.size() == maxs.size());
        DataKey = std::make_pair(s, maxs);
        nKeys = s.size();
        Meta = std::vector<std::string>(s.size(), "null");
        CurrentKey = "";
    }

    // @brief Include @meta description of the @s
    DataKeys(const std::vector<std::string>& s, const std::vector<uint>& maxs, const std::vector<std::string>& meta) {
        assert(s.size() == maxs.size());
        assert(s.size() == meta.size());
        DataKey = std::make_pair(s, maxs);
        nKeys = s.size();
        Meta = meta;
        CurrentKey = "";
    }

    std::string MakeKey(const std::vector<int> &idx) {
        assert(idx.size() == DataKey.first.size());
        std::string key = "";
        for (size_t i = 0; i  < DataKey.first.size(); ++i) 
           key += (DataKey.first[i] + std::to_string(idx[i])) + "-";
        key += "-";
        CurrentKey = key;
        return key;
    }

    void Print() {
       std::cout << "\n----- Data key summary -----" << std::endl;
       for (size_t i = 0; i < nKeys; ++i) {
          if (DataKey.second[i] == 1)
              std::cout << DataKey.first[i] << "[1] : " << Meta[i] << "." << std::endl;
          else 
              std::cout << DataKey.first[i] << "[1-" << DataKey.second[i] << "] : " << Meta[i] << "." << std::endl;
       }
       std::cout << "Current key = " << CurrentKey << std::endl;
       std::cout << "----------------------------\n" << std::endl;
    }

    std::pair<std::vector<std::string>, std::vector<uint>> DataKey;
    std::string CurrentKey;
    std::vector<std::string> Meta;
    size_t nKeys;
};


int main(int argc, char *argv[])
{
  // ======================= BOILER PLATE START ==========================
  std::cout << std::setprecision(std::numeric_limits<float_>::digits10) << std::showpoint;
  std::vector<int> KeyIndices;
  for (int i = 2; i < argc; ++i) 
      KeyIndices.push_back(std::atoi(argv[i]));
  //DataKeys DK({"s", "b", "p"}, {3, 9, 1}, {"System", "Geometry parameter: Internal angle", "Partitions" });
  DataKeys DK({"b", "p"}, {2, 2}, {"Geometry parameter: Internal angle", "Partitions" });
  std::string Key = DK.MakeKey(KeyIndices);
  DK.Print();
  // ======================= BOILER PLATE END ==========================

  // ======================= BILLIARD DEFINITION START ==========================
  std::vector<std::unique_ptr<Obstacle<float_>>> Obstacles(1);

  // Define rational internal angles
  //std::vector<float_> alpha_rational = seq<float_>((boost::math::constants::pi<float_>() / (float_)4), (3 * M_PI / (float_)4), DK.DataKey.second[0]);
  const float_ GR = boost::math::constants::phi<float_>(); // Golden ratio
  //std::vector<float_> alpha_rational   = seq<float_>(boost::math::constants::pi<float_>()_2, (boost::constants::pi<float_>()/((float_)2)) / (float_)5, DK.DataKey.second[0]);
  //std::vector<float_> alpha_irrational = seq<float_>((5 * boost::math::constants::pi<float_>()) / (2 + 5 * GR), (5 * M_PI) / (25 * GR + 10), DK.DataKey.second[0]);


  std::vector<float_> alpha_rational   = seq<float_>(boost::math::constants::pi<float_>()_2, (boost::constants::pi<float_>()/((float_)2)) / (float_)5, 10);
  std::vector<float_> alpha_irrational = seq<float_>((5 * boost::math::constants::pi<float_>()) / (2 + 5 * GR), (5 * M_PI) / (25 * GR + 10), 10);

  std::vector<float_> alpha = { alpha_rational[2], alpha_irrational[2] };

  assert((size_t)KeyIndices[0] < alpha_rational.size() - 1);
  assert((size_t)KeyIndices[1] <= DK.DataKey.second[1]);

  Polygon<float_> Channel = getParallelChannel<float_>(alpha[KeyIndices[0] - 1], (1 / (float_)2), -(1 / (float_)10));
  Channel.Write("Data/BilliardTable.dat");
  Obstacles[0] = std::make_unique<Polygon<float_>>(Channel);
  Billiard<float_> Table(Obstacles);
  // ======================= BILLIARD DEFINITION END ==========================


  Simulation<float_> Sim;

  // As many particles as possible to escape within the time limit
  Sim = Simulation<float_> {
      .Time = getExponentialTimes<float_>((float_)1e6, 1000),
      .L = seq<float_>((float_)10.0, (float_)50.0, 5),
      .Symbols = {"."},
      .DataDestination = "Data/s0-" + Key,
      .RunTime = ".",
      .MomentStep = 1,
      .TotalTime = (float_)1e6,
      .nInitialConditions = 100000,
      .nCollisions = 0,
      .nMoments = 4};
  Sim.Write();
  EscapeToAbsorbingBoundary(Table, Sim);
  //std::thread th1(EscapeToAbsorbingBoundary<float_>, std::ref(Table), std::ref(Sim));

  //// Ensures all particles escape to L
  for (size_t i = 0; i < 4; ++i) {
    float_ TT = std::pow(10, i+3);
    Sim = Simulation<float_> {
        .Time = getExponentialTimes<float_>(TT, 1000),
        .L = seq<float_>((float_)10.0, (float_)50.0, 5),
        .Symbols = {"."},
        .DataDestination = "Data/s1-t" + std::to_string(i+1) + "-" + Key,
        .RunTime = ".",
        .MomentStep = 1,
        .TotalTime = TT,
        .nInitialConditions = 100000,
        .nCollisions = 0,
        .nMoments = 4};
    ForcedEscapeToAbsorbingBoundary(Table, Sim);
  }
  //std::thread th2(ForcedEscapeToAbsorbingBoundary<float_>, std::ref(Table), std::ref(Sim));
  
  Sim = Simulation<float_> {
      .Time = getExponentialTimes<float_>((float_)1e6, 1000),
      .L = {0},
      .Symbols = {"0404040404", "1313131313", "0431043104", "1340134013"},
      .DataDestination = "Data/s2-" + Key,
      .RunTime = ".",
      .MomentStep = 1,
      .TotalTime = (float_)1e5,
      .nInitialConditions = 10000,
      .nCollisions = 10000,
      .nMoments = 4};
  EscapeToSymbolicSequence(Table, Sim);
  //std::thread th3(EscapeToSymbolicSequence<float_>, std::ref(Table), std::ref(Sim));


  //th1.join();
  //th2.join();
  //th3.join();



  return 0;
}




    // New rational parameters
    // const vector<int> dy_map = {0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 1, 2, 3, 4, 5, 6, 7, 8};
    // // const vector<int> d_map = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    // const vector<int> d_map = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    // FLOAT Delta_x = 1 / 2.0L;
    // // vector<FLOAT> alpha = seq(boost::math::constants::pi<float_>()_2, (boost::constants::pi<float_>()/((float_)2)) / 5.0L, 27);
    // FLOAT eps = boost::math::constants::pi<float_>() / 180;
    // // vector<FLOAT> alpha = seq(boost::math::constants::pi<float_>()_2 - eps, (boost::constants::pi<float_>()/((float_)2)) + eps, 5);
    // vector<FLOAT> alpha = {boost::math::constants::pi<float_>()_2};
    // reverse(alpha.begin(), alpha.end());
    // vector<FLOAT> dy(alpha.size());
    // for (size_t i = 0; i < alpha.size(); ++i)
    // {
    //     FLOAT val = Delta_x * cosl(alpha[i] / 2.0L) / sinl(alpha[i] / 2.0L);
    //     std::cout << alpha[i] * 180 / boost::math::constants::pi<float_>() << " " << val << std::endl;
    //     dy[i] = val;
    // }

    // // FLOAT Delta_y = dy[dy_map[BID]];
    // FLOAT Delta_y = dy[BID];
    // FLOAT d = Delta_y - Delta_y / 2.0L;

    // FLOAT d = -1;
    // if (d_map[BID] == 0)
    //     d = Delta_y + Delta_y / 2.0L;
    // else if (d_map[BID] == 1)
    //     d = Delta_y - Delta_y / 2.0L;

    // New irrational parameters
    // const vector<int> dy_map = {0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 1, 2, 3, 4, 5, 6, 7, 8};
    // const vector<int> d_map = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    // FLOAT Delta_x = 1 / 2.0L;
    // const FLOAT GR = (sqrtl(5) + 1.0L) / 2.0L;
    // vector<FLOAT> alpha = seq((5.0L * boost::math::constants::pi<float_>()) / (2.0L + 5.0L * GR), (5.0L * M_PI) / (25.0L * GR + 10.0L), 9);
    // reverse(alpha.begin(), alpha.end());
    // vector<FLOAT> dy(alpha.size());
    // for (size_t i = 0; i < alpha.size(); ++i)
    // {
    //     FLOAT val = Delta_x * cosl(alpha[i] / 2.0L) / sinl(alpha[i] / 2.0L);
    //     std::cout << alpha[i] << " " << val << std::endl;
    //     dy[i] = val;
    // }

    // FLOAT Delta_y = dy[dy_map[BID]];
    // FLOAT d = -1;
    // if (d_map[BID] == 0)
    //     d = Delta_y + Delta_y / 2.0L;
    // else if (d_map[BID] == 1)
    //     d = Delta_y - Delta_y / 2.0L;



    // Symbolic escape
    //Trajectory<float_> DiscreTraj = Table.getTrajectory(DC, N);
    //std::vector<float_> FPTSymb = DiscreTraj.getEscapeTimes(symb);

    // Particle<float_> IC = Table.getInitialState(0, -1, T);
    // IC.Pos = Point<float_>(0, (Channel.getVertices()[1].y + Channel.getVertices()[4].y) / float_(2));
    // IC.Vel = Point<float_>(cos(theta[i]), sin(theta[i]));
    // IC.WriteState(DataDir + "InitialConditions.dat");
