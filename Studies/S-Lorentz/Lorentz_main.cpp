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
  std::cout << "Arugments received: ";
  for (int i = 0; i < argc; ++i)
    std::cout << argv[i] << " ";
  std::cout << std::endl;

  // ======================= BOILER PLATE START ==========================
  std::cout << std::setprecision(std::numeric_limits<float_>::digits10) << std::showpoint;
  std::vector<int> KeyIndices;
  for (int i = 2; i < argc; ++i)
    KeyIndices.push_back(std::atoi(argv[i]));
  DataKeys DK({"p", "b", "d"}, {100, 6, 2}, {"Partitions", "Boundary parameter: [0-2] Irrational, [3-5] Rational", "Pore widths: [0] Finite horizon, [1] Infinite horizon"});
  std::string Key = DK.MakeKey(KeyIndices);
  DK.Print();
  // ======================= BOILER PLATE END ==========================

  // ======================= BILLIARD DEFINITION START ==========================
  std::vector<std::unique_ptr<Obstacle<float_>>> Obstacles(2);

  // Fundamental cell
  Polygon<float_> Square({(float_)0, (float_)1, (float_)1, (float_)0}, {(float_)0, (float_)0, (float_)1, (float_)1}, {1, 3, 0, 2}, true);

  float_ radius = 2 / ((float_)10);
  Circle<float_> Scatterer(radius, Point<float_>((float_)0.5, (float_)0.5));

  Obstacles[0] = std::make_unique<Polygon<float_>>(Square);
  Obstacles[1] = std::make_unique<Circle<float_>>(Scatterer);

  Billiard<float_> Table(Obstacles);
  // ======================= BILLIARD DEFINITION END ==========================
  //
  float_ T = (float_)1e5;
  std::vector<float_> Time = getExponentialTimes<float_>(T, 1000);

  std::vector<float_> stime = {(float_)100, (float_)1000, (float_)10000, (float_)100000};
  std::vector<float_> TIME = Time;
  TIME.insert(TIME.end(), stime.begin(), stime.end());
  sort(TIME.begin(), TIME.end());

  WriteVector<float_>("Data/TIME.dat", TIME);

  Simulation<float_> Sim{
      .Time = TIME,
      .L = {0},
      .Symbols = {""},
      .DataDestination = "Data/" + Key,
      .RunTime = ".",
      .MomentStep = (float_)0.5,
      .TotalTime = TIME.back(),
      .nInitialConditions = 20000,
      .nCollisions = 0,
      .nMoments = 20};
  //Sim.Write();

  //WriteTrajectory_All<float_>(Table, Sim);
  //VelocityCorrelations<float_>(Table, Sim);
  ContinuousSystemSummary2D<float_>(Table, Sim);

  return 0;
}

// #include <iostream>
// #include <algorithm>
// #include <chrono>
// #include <thread>
// #include <numeric>
// #include "Billiard.h"
// #include "Trajectory.h"
// #include "OnlineStatistics.h"
// #include "Misc.h"
// #include "Ensemble.h"
// #include "Simulation.h"

// typedef float_ldub float_;

// const int nDigits = std::numeric_limits<float_>::digits10;

// int main(int argc, char *argv[])
// {

//   // ======================= BOILER PLATE START ==========================
//   if (argc != 2 && argc != 3)
//   {
//     std::cout << "Not enough command line arguments." << std::endl;
//     std::cout << "Require:    Job ID" << std::endl;
//     exit(EXIT_FAILURE);
//   }
//   std::string DataFilePrefix;
//   if (argc == 2)
//     DataFilePrefix = std::to_string(std::atoi(argv[1])) + "-";
//   else
//     DataFilePrefix = std::to_string(std::atoi(argv[1])) + "-" + std::to_string(std::atoi(argv[2])) + "-";

//   const int JobID = std::atoi(argv[1]);
//   std::cout << "Running job " << argv[1] << " on partition " << (argc == 2 ? "0" : argv[2]) << "." << std::endl;
//   std::string DataDir = "DATA/" + DataFilePrefix;
//   std::cout << "Saving data to: " << DataDir << std::endl;
//   std::cout << "No. precision digits: " << std::numeric_limits<float_>::digits10 << std::endl;
//   std::cout << std::setprecision(std::numeric_limits<float_>::digits10) << std::showpoint;
//   // ======================= BOILER PLATE END ==========================

//   // ======================= BILLIARD DEFINITION START ==========================
//   std::vector<std::unique_ptr<Obstacle<float_>>> Obstacles(1);

//   float_ PHI = boost::math::constants::phi<float_>();
//   float_ PI = boost::math::constants::pi<float_>();

//   std::vector<float_> alpha_irrational = seq((5 * PI) / (2 + 5 * PHI), (5 * PI) / (25 * PHI + 10), 3); // IRRATIONAL
//   // std::vector<float_> alpha_rational = seq(PI / (float_)2, PI / (float_)10, 9);                        // RATIONAL
//   std::vector<float_> alpha_rational = seq(PI / (float_)2, PI / (float_)10, 9); // RATIONAL

//   int jid = JobID;
//   float_ alpha;
//   if ((size_t)jid > alpha_irrational.size() - 1)
//   {
//     jid = jid % alpha_rational.size();
//     alpha = alpha_rational[jid];
//   }
//   else // irrationals 0-8
//   {
//     alpha = alpha_irrational[jid];
//   }
//   std::vector<float_> pw = pseq<float_>(-0.5, 0, 5);
//   // Polygon<float_> Channel = getParallelChannel<float_>(alpha, 1 / static_cast<float_>(2), pw[jid]);
//   Polygon<float_> Channel = getParallelChannel<float_>(PI / 2, 1 / static_cast<float_>(2), pw[jid]);
//   // Polygon<float_> Channel = getParallelChannel<float_>(alpha, 1 / ((float_)(2)), (float_)0);
//   Channel.Write(DataDir + "Channel.dat");
//   Obstacles[0] = std::make_unique<Polygon<float_>>(Channel);
//   Billiard<float_> Table(Obstacles);
//   // ======================= BILLIARD DEFINITION END ==========================

//   // ======================= SIMULATION DETAILS ==========================

//   float_ T = 1e2; // Total time / max time
//   int M = 10000;  // Initial conditions
//   int N = 10;     // Number of collisions
//   int nMoments = 5;
//   float_ MomentStep = static_cast<float_>(1);
//   std::vector<float_> Time = getExponentialTimes<float_>(T, 100);
//   // std::vector<float_> Time = getLinearTimes<float_>(T, 100);
//   // std::vector<float_> Time = getLinearTimes<float_>(T, 1);
//   std::vector<float_> L = {10};

//   std::vector<std::string> symb = {"040404040", "04040404040", "0404040404040", "040404040404040"};

//   Simulation<float_> Sim{
//       .Time = Time,
//       .L = L,
//       .Symbols = symb,
//       .DataDestination = DataDir,
//       .RunTime = "",
//       .MomentStep = MomentStep,
//       .TotalTime = Time.back(),
//       .nInitialConditions = M,
//       .nCollisions = N,
//       .nMoments = nMoments};

//   // PrecisionTest<float_>(Table, Sim);

//   // PrecisionTest<float_50>(Table, Sim);

//   // ContinuousSystemSummary<float_>(Table, Sim);
//   // VelocityCorrelations<float_>(Table, Sim);
//   // ForcedEscapeToAbsorbingBoundary<float_>(Table, Sim);
//   // FirstPassageTimes<float_>(Table, Sim);
//   EscapeToAbsorbingBoundary<float_>(Table, Sim);
//   // MeanFirstPassage<float_>(Table, Sim);
//   // // CollisionSystemSummary<float_>(Table, Sim);
//   // // WriteSampleTrajectory(Table, Sim);
//   exit(EXIT_FAILURE);

//   // if (L[0] < 0.5)
//   //   std::cout << "Warning: Escape times will be checked on fundamental domain." << std::endl;
//   // WriteVector(DataDir + "L.dat", L);

//   // std::vector<float_> theta = seq<float_>(0.0, 2 * boost::math::constants::pi<float_>(), M);

//   // std::vector<float_> PosSum(L.size(), float_(0));
//   // std::vector<float_> AbsPosSum(L.size(), float_(0));
//   // std::vector<Ensemble<float_>> MeanXPosSum(L.size(), Ensemble<float_>(N)); // Over partial sums

//   // Ensemble<float_> MeanEscapeTime(L.size());
//   // // std::vector<Ensemble<float_>> FiniteMoments(2, Ensemble<float_>(Time.size()));

//   // std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
//   // for (int i = 0; i < M; ++i)
//   // {
//   //   if (10 * i % M == 0 && i != 0)
//   //     std::cout << "Simulated " << i << " initial conditions." << std::endl;

//   //   // Particle<float_> IC = Table.getInitialState(0, -1, T);
//   //   // IC.Pos = Point<float_>(0, (Channel.getVertices()[1].y + Channel.getVertices()[4].y) / float_(2));
//   //   // IC.Vel = Point<float_>(cos(theta[i]), sin(theta[i]));
//   //   // IC.WriteState(DataDir + "InitialConditions.dat");

//   //   // For finite trajectories we need discrete trajectories
//   //   EscapeData<float_> ED(L);
//   //   Particle<float_> IC;
//   //   // while (ED.EscapeTimes.back() == 0)
//   //   // {
//   //   IC = Table.getInitialState(0, -1, T);
//   //   Trajectory<float_> FiniteTraj = Table.getFiniteTrajectory(IC, N, L.back());
//   //   // Trajectory<float_> FiniteTraj = Table.getFiniteTrajectory(IC, L.back());
//   //   if (Table.Failed)
//   //   {
//   //     Table.ErrorReport(IC, i);
//   //     i--; // Want M initial conditions
//   //   }
//   //   ED = FiniteTraj.getEscapeData(L);
//   //   // }

//   //   for (size_t k = 0; k < L.size(); ++k)
//   //   {
//   //     std::vector<float_> pos = FiniteTraj.getPositionsUpToEscape(L[k]);
//   //     // std::vector<float_> pos = FiniteTraj.getDiffsUpToEscape(L[k]);

//   //     PosSum[k] = std::accumulate(pos.begin(), pos.end(), (float_)0);

//   //     for (auto &f : pos)
//   //       f = f < 0 ? -f : f;

//   //     std::vector<float_> PartialSum(pos.size());
//   //     std::partial_sum(pos.begin(), pos.end(), PartialSum.begin());
//   //     MeanXPosSum[k].Push(PartialSum);

//   //     AbsPosSum[k] = std::accumulate(pos.begin(), pos.end(), (float_)0);
//   //   }

//   //   // Finite moments
//   //   // Trajectory<float_> FiniteContinuousTraj = FiniteTraj;
//   //   // FiniteContinuousTraj.MakeContinuous(Time);
//   //   // std::vector<std::vector<float_>> DXp = FiniteContinuousTraj.getDisplacements_x(1, 2);
//   //   // for (int j = 0; j < 2; ++j)
//   //   //   FiniteMoments[j].Push(DXp[j]);

//   //   // WriteVector(DataDir + "x.dat", getCoordFromVectorPoint(FiniteTraj.getPositions(), "x"));

//   //   MeanEscapeTime.Push(ED.EscapeTimes);
//   //   // Symbolic escape
//   //   // Trajectory<float_> DiscreTraj = Table.getTrajectory(DC, N);
//   //   // std::vector<float_> FPTSymb = DiscreTraj.getEscapeTimes(symb);

//   //   // Finite moments
//   //   // Trajectory<float_> FiniteContinuousTraj = FiniteTraj;
//   //   // FiniteContinuousTraj.MakeContinuous(Time);
//   //   // fMSD.Push(FiniteContinuousTraj.getDisplacements_x2());
//   //   WriteVector(DataDir + "xPosSum.dat", PosSum);
//   //   WriteVector(DataDir + "AbsxPosSum.dat", AbsPosSum);
//   //   WriteVector(DataDir + "EscapeTimes.dat", ED.EscapeTimes);
//   //   WriteVector(DataDir + "EscapeAbscissa.dat", ED.EscapePosition_x);
//   //   WriteVector(DataDir + "EscapeDirection.dat", ED.EscapeDirection);
//   //   IC.WriteState(DataDir + "InitialConditions.dat");
//   // }
//   // WriteVector(DataDir + "MeanEscapeTimes.dat", MeanEscapeTime.PullMeans());

//   // // for (size_t i = 0; i < FiniteMoments.size(); ++i)
//   // //   FiniteMoments[i].Write(Sim.DataDestination + "DS_FiniteMoment_" + std::to_string(i) + ".dat");

//   // for (size_t k = 0; k < L.size(); ++k)
//   //   WriteVector(DataDir + "MeanXPosSum.dat", MeanXPosSum[k].PullMeans());

//   // std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
//   // const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";
//   // Sim.Write();
//   // Sim.Print(RunTime);

//   return 0;
// }
