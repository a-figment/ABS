#include <iostream>
#include <algorithm>
#include <chrono>
#include "Billiard.h"
#include "Trajectory.h"
#include "OnlineStatistics.h"

typedef float_ldub float_;

struct SimulationDetails
{
  int nInitialConditions;
  float_ nTimeUnits;
  std::string runtime;

  void Print()
  {
    std::cout << nInitialConditions << " initial conditions." << std::endl;
    std::cout << nTimeUnits << " time units." << std::endl;
    std::cout << runtime << " time units." << std::endl;
  }
};

const int nDigits = std::numeric_limits<float_>::digits10;

template <typename T>
std::vector<T> getExponentialTimes(const T &totalTime, const size_t &nSteps)
{
  std::cout << "Using log time with " << nSteps << " steps." << std::endl;
  T step = log(totalTime);
  std::vector<T> Time(nSteps);
  for (size_t i = 0; i < nSteps; ++i)
  {
    T val = (step / (nSteps - static_cast<T>(1))) * static_cast<T>(i);
    Time[i] = exp(val);
  }
  return Time;
}

template <typename T>
std::vector<T> getLinearTimes(const T &totalTime, const T &dt)
{
  std::cout << "Using linear times: ." << std::endl;
  std::vector<T> Time(static_cast<size_t>(totalTime)); // Continuous times to observe particle
  for (size_t i = 1; i <= Time.size(); ++i)
    Time[i - 1] = dt * static_cast<T>(i);
  return Time;
}

int main(int argc, char *argv[])
{
  if (argc != 2)
  {
    std::cout << "Not enough command line arguments." << std::endl;
    exit(EXIT_FAILURE);
  }
  std::cout << "Running job " << argv[1] << std::endl;
  std::string DataDir = "../Studies/S-MPTest/Data/" + std::to_string(std::atoi(argv[1])) + "-";
  std::cout << "Saving data to: " << DataDir << std::endl;
  std::cout << std::setprecision(std::numeric_limits<float_>::max_digits10) << std::showpoint;

  uint nObstacles = 1;
  std::vector<std::unique_ptr<Obstacle<float_>>> Obstacles(nObstacles);

  const float_ Delta_x = 1 / static_cast<float_>(2);
  float_ GR = boost::math::constants::phi<float_>();
  std::vector<float_> S_i, DY, ALPHA;
  for (int i = 3; i >= 0; i--)
  {
    float_ si = (13 / static_cast<float_>(126)) * i + (6 / static_cast<float_>(7));
    float_ alpha = (si * boost::math::constants::pi<float_>()) / GR;
    float_ delta_y = Delta_x * (cos(alpha / static_cast<float_>(2)) / sin(alpha / static_cast<float_>(2)));
    // std::cout << std::setprecision(std::numeric_limits<float_>::max_digits10) << std::showpoint;
    // std::cout << "(s_" << i << ", ALPHA, DY) = (" << si << ", " << alpha << ", " << delta_y << ")" << std::endl;
    S_i.push_back(si);
    ALPHA.push_back(alpha);
    DY.push_back(delta_y);
  }

  float_ dy = DY[3];
  // float_ dy = 1 / static_cast<float_>(2);
  float_ d = dy; // / static_cast<float_>(1.5);
  std::cout << std::setprecision(std::numeric_limits<float_>::digits10) << std::showpoint;
  std::cout << "(d, dy) = " << Point<float_>(d, dy) << std::endl;

  // Define channel
  Polygon<float_> Channel({-Delta_x, float_(0), Delta_x, Delta_x, float_(0), -Delta_x}, {float_(0), dy, float_(0), d, d + dy, d}, {2, 5}, false);
  Channel.Write(DataDir + "Channel.dat");
  Obstacles[0] = std::make_unique<Polygon<float_>>(Channel);
  Billiard<float_> Table(Obstacles);

  // ===============
  float_ T = 1e6; // Total time
  int M = 1;      // Initial conditions

  std::vector<float_> Time = getExponentialTimes<float_>(T, 1000);
  std::vector<OnlineStatistics<float_>> MSD(Time.size(), OnlineStatistics<float_>());

  std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < M; ++i)
  {
    if (10 * i % M == 0 && i != 0)
      std::cout << "Simulated " << i << " initial conditions." << std::endl;

    Particle<float_> IC = Table.getInitialState(0, -1, T);
    IC.WriteState(DataDir + "InitialConditions.dat");
    Particle<float_> Pt = IC;

    std::chrono::high_resolution_clock::time_point z0 = std::chrono::high_resolution_clock::now();
    Trajectory<float_> ContTraj = Table.getContinuousTrajectory(Pt, Time);
    std::chrono::high_resolution_clock::time_point z1 = std::chrono::high_resolution_clock::now();
    std::cout << "\033[1;32mTrajectory run time: " << std::chrono::duration_cast<std::chrono::hours>(z1 - z0).count() << " hours = " << std::chrono::duration_cast<std::chrono::minutes>(z1 - z0).count() << " minutes = " << std::chrono::duration_cast<std::chrono::milliseconds>(z1 - z0).count() << " milliseconds. \n\033[0m" << std::endl;

    if (Table.Failed)
    {
      Table.ErrorReport(Pt, i);
      i--; // Want M initial conditions
    }
    std::vector<Point<float_>> p = ContTraj.getPositions();
    for (size_t j = 0; j < p.size(); ++j)
    {
      float_ dx = fabs(p[j].x - p[0].x);
      MSD[j].Update(dx * dx);
    }
  }

  std::ofstream X(DataDir + "MSD-" + std::to_string(nDigits) + "-digits.dat", std::ios::out);
  X << std::setprecision(std::numeric_limits<float_>::digits10) << std::showpoint;
  for (size_t i = 0; i < MSD.size(); ++i)
    X << MSD[i].getMean() << " ";
  X << "\n";

  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  std::cout << "\033[1;32mSimulation run time: " << std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count() << " hours = " << std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count() << " minutes = " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << " milliseconds. \n\033[0m" << std::endl;

  return 0;
}
