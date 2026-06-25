#ifndef ESCAPEDATA_H_INCLUDED
#define ESCAPEDATA_H_INCLUDED

#include "Trajectory.h"
#include <sstream>

// @brief Upon an escape event, various states of the particle can be stored.
template <typename T>
struct EscapeData
{
    EscapeData(const std::vector<T> &EscapeEvents, const Trajectory<T> &traj)
    {
        Traj = traj;
        EscapeTimes.resize(EscapeEvents.size());
        EscapePosition_x.resize(EscapeEvents.size());
        EscapePosition_y.resize(EscapeEvents.size());
        EscapeDirection.resize(EscapeEvents.size());
        xPositionSum.resize(EscapeEvents.size());
        nEscapers.resize(EscapeEvents.size());
    }

    void setEscapeData(const std::vector<T> &L);
    void setEscapeData(const std::vector<std::string> &Symbols);

    std::vector<T> EscapeTimes;
    std::vector<T> EscapePosition_x;
    std::vector<T> EscapePosition_y;
    std::vector<T> EscapeDirection;
    std::vector<T> xPositionSum; // Sum of x-positions up until an escape
    std::vector<int> nEscapers;  // Sum of x-positions up until an escape

    Trajectory<T> Traj;
};

template <typename T>
void EscapeData<T>::setEscapeData(const std::vector<T> &L)
{

    std::vector<Point<T>> Positions = Traj.getLiftedPositions();
    std::vector<Point<T>> Velocities = Traj.getVelocities();
    std::vector<T> Time = Traj.getTime();
    size_t startidx = 0;

    for (size_t i = 0; i < L.size(); ++i)
    {
        for (size_t j = startidx; j < Positions.size(); ++j)
        {
            if (fabs(Positions[j].x) > L[i]) // Escape occurred
            {
                T l = L[i];
                if (Velocities[j - 1].x < static_cast<T>(0))
                    l = -L[i];
                Point<T> EPos = getPointOnVerticalLine(Positions[j - 1], Velocities[j - 1].y / Velocities[j - 1].x, l);
                const T d = Distance2D<T>(Positions[j - 1], EPos);

                // Escape data
                EscapeTimes[i] = Time[j - 1] + d;
                EscapePosition_x[i] = EPos.x;
                EscapePosition_y[i] = EPos.y;
                EscapeDirection[i] = atan2(Velocities[j - 1].y, Velocities[j - 1].x);

                nEscapers[i]++;

                startidx = j - 1;
                break;
            }
        }
    }
}

template <typename T>
void EscapeData<T>::setEscapeData(const std::vector<std::string> &Symbols)
{
    std::vector<int> IDs = Traj.getSegmentIDs();

    // Convert sequence to string
    std::ostringstream oss;
    std::copy(IDs.begin() + 1, IDs.end(), std::ostream_iterator<int>(oss, ""));
    std::string sID = oss.str();

    for (size_t i = 0; i < Symbols.size(); ++i)
    {
        size_t idx = sID.find(Symbols[i], 0);
        if (idx != std::string::npos)
        {
            EscapeTimes[i] = Traj.MakeParticle(idx).TotalTime;
            EscapePosition_x[i] = Traj.MakeParticle(idx).Pos.x;
            EscapePosition_y[i] = Traj.MakeParticle(idx).Pos.y;
            EscapeDirection[i] = atan2(Traj.MakeParticle(idx).Vel.y, Traj.MakeParticle(idx).Vel.x);
        }
        else
        {
            nEscapers[i]++;
        }
    }
}

// template <typename T>
// std::vector<T> EscapeData<T>::getEscapeTimesAtRandomCollision(const int &nToGet)
// {
//     std::vector<T> ETARC(nToGet, -1);

//     std::random_device rd;  //Will be used to obtain a seed for the random number engine
//     std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
//     std::uniform_int_distribution<> distrib(1, static_cast<int>(Particles.size()));

//     Random<int> RInt2(1, static_cast<int>(Particles.size()));
//     for (int i = 0; i < nToGet; ++i)
//     {
//         int RandomCollision = distrib(gen);
//         ETARC[i] = Particles[RandomCollision].TotalTime;
//     }
//     return ETARC;
                // std::vector<T>::iterator upTo{Positions.begin() + j - 1};

                // xPositionSum[i] = std::accumulate(Positions.begin(), Positions.begin() + (j - 1), (T)0) + L[i];
                // for (auto x : Positions)
                // xPositionSum[i] += abs(x.x);
                // std::transform_reduce(Positions.begin(), vec.end(),0,std::plus<>{}, static_cast<double (*)(double)>(std::fabs));
// }

#endif
