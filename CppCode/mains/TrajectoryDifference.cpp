#include <iostream>
#include <algorithm>
#include <chrono>
#include "Billiard.h"
#include "Trajectory.h"
#include "OnlineStatistics.h"

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

template <typename T>
Trajectory<T> getTraj()
{
    uint nObstacles = 1;

    const T Delta_x = 1 / static_cast<T>(2);
    const T alpha = ((13 / static_cast<T>(126)) * T(3) + (6 / static_cast<T>(7)) * boost::math::constants::pi<T>()) / (boost::math::constants::phi<T>());
    const T dy = Delta_x * (cos(alpha / static_cast<T>(2)) / sin(alpha / static_cast<T>(2)));
    const T d = dy;
    Polygon<T> Channel({-Delta_x, T(0), Delta_x, Delta_x, T(0), -Delta_x}, {T(0), dy, T(0), d, d + dy, d}, {2, 5}, false);

    std::vector<std::unique_ptr<Obstacle<T>>> Obstacles(nObstacles);
    Obstacles[0] = std::make_unique<Polygon<T>>(Channel);
    Billiard<T> Table_50(Obstacles);

    std::cout << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;
    std::cout << "(d, dy) = " << Point<T>(d, dy) << std::endl;

    std::vector<T> x = {T(0)};
    std::vector<T> y = {dy + T(0.1)};

    Point<T> eps(T(10) * std::numeric_limits<T>::epsilon(), T(10) * std::numeric_limits<T>::epsilon());
    Point<T> ipos(x[0], y[0]);
    Point<T> vel = Channel.getVertices()[1] - ipos;
    Particle<T> Pt(ipos, vel / Norm(vel) + eps, -1, INFINITY);

    Trajectory<T> DiscTraj = Table_50.getDiscreteTrajectory(Pt, 100000);
    return DiscTraj;
}

int main()
{

    Trajectory<float_50> T50 = getTraj<float_50>();
    Trajectory<float_200> T200 = getTraj<float_200>();

    // Point<float_50> p50 = T50.getPositions();
    // Point<float_200> p200 = T200.getPositions();
    auto p50 = T50.getPositions();
    auto p200 = T200.getPositions();

    std::vector<float_200> diffs(T50.getPositions().size());
    for (size_t i = 0; i < T50.getPositions().size(); ++i)
    {
        diffs[i] = fabs(p50[i].x - p200[i].x);
    }

    std::ofstream X("diffs.dat", std::ios::out);
    X << std::setprecision(std::numeric_limits<float_200>::digits10) << std::showpoint;
    for (size_t i = 0; i < diffs.size(); ++i)
        X << diffs[i] << " ";
    X << "\n";
    // DiscTraj.Write(DataDir + "dTraj-" + std::to_string(std::numeric_limits<float_800>::digits10) + "-digits.dat");

    return 0;
}
