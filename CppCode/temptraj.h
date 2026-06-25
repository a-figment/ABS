#ifndef TRAJECTORY_H_INCLUDED
#define TRAJECTORY_H_INCLUDED

#include "Particle.h"
#include "Random.h"
#include "Ensemble.h"
#include <sstream>

// @brief A trajectory is a collection of particle properties. We take a data-oriented approach and store
//        each property as a simple vector whose size is determined at run time. Optionally, if the size isn't known
//        one can just as easily reserve a rough estimate and use @Push_ instead of @Push.
// @Push() - Unpacks a particle object into vectors of each fundamental type (FT). We consider FT's as
//                  1. Position (x,y)
//                  2. Velocity (x,y)
//                  3. Segmend identifiers
//                  4. Obstacle identifiers
//
template <typename T>
class Trajectory
{
public:
    Trajectory() = default;

    // @Length - Total number of particles to consider
    Trajectory(const size_t &Size)
    {
        nCurrentElements = 0;
        nTotalElements = Size;
    }

    void Reserve(const size_t &EstimatedSize)
    {
        Positions.reserve(EstimatedSize);
        Velocities.reserve(EstimatedSize);
        FreePaths.reserve(EstimatedSize);
        SegmentIDs.reserve(EstimatedSize);
        ObstacleIDs.reserve(EstimatedSize);
    }

    Trajectory(const size_t &Size, const Particle<T> &IP)
    {
        nCurrentElements = 0;
        nTotalElements = Size;

        Positions = std::vector<Point<T>>(Size, Point<T>());
        CellIDs = std::vector<Point<cid_int>>(Size, Point<cid_int>());
        Velocities = std::vector<Point<T>>(Size, Point<T>());
        FreePaths = std::vector<T>(Size, (T)0);
        SegmentIDs = std::vector<int>(Size, -2); // -1 => particle is floating. -2 => it is undefined
        ObstacleIDs = std::vector<int>(Size, -2);

        Push(IP);
    }

    void Push(const Particle<T> &P)
    {
        Positions[nCurrentElements] = P.Pos; //+ Point<T>((T)P.CellID.x, (T)P.CellID.y); // Lifted position
        CellIDs[nCurrentElements] = P.CellID;
        Velocities[nCurrentElements] = P.Vel;
        FreePaths[nCurrentElements] = P.FreePathLength;
        SegmentIDs[nCurrentElements] = P.SegmentID;
        ObstacleIDs[nCurrentElements] = P.ObstacleID;
        nCurrentElements++;
    }

    void Push_(const Particle<T> &P)
    {
        Positions.push_back(P.Pos); // + Point<T>((T)P.CellID.x, (T)P.CellID.y)); // Lifted position
        CellIDs.push_back(P.CellID);
        Velocities.push_back(P.Vel);
        FreePaths.push_back(P.FreePathLength);
        SegmentIDs.push_back(P.SegmentID);
        ObstacleIDs.push_back(P.ObstacleID);
        nCurrentElements++;
    }

    void Pop(const size_t &offset)
    {
        Positions.erase(Positions.begin() + offset, Positions.end());
        CellIDs.erase(CellIDs.begin() + offset, CellIDs.end());
        Velocities.erase(Velocities.begin() + offset, Velocities.end());
        FreePaths.erase(FreePaths.begin() + offset, FreePaths.end());
        SegmentIDs.erase(SegmentIDs.begin() + offset, SegmentIDs.end());
        ObstacleIDs.erase(ObstacleIDs.begin() + offset, ObstacleIDs.end());
        nTotalElements = Positions.size();
    }

    // Base types
    size_t getNumElements() const { return nCurrentElements; }
    size_t getLength() const { return nTotalElements; }

    std::vector<Point<T>> getPositions() const { return Positions; }
    std::vector<Point<T>> getLiftedPositions() const
    {
        std::vector<Point<T>> LiftedPositions(Positions.size(), Point<T>());
        for (size_t i = 0; i < Positions.size(); ++i)
            LiftedPositions[i] = Point<T>(Positions[i].x + (T)CellIDs[i].x, Positions[i].y + (T)CellIDs[i].y);
        return LiftedPositions;
    }
    std::vector<Point<cid_int>> getCellIDs() const { return CellIDs; }
    std::vector<Point<T>> getVelocities() const { return Velocities; }
    std::vector<T> getFreePaths() const { return FreePaths; }
    std::vector<int> getSegmentIDs() const { return SegmentIDs; }

    // Coordinate access to point objects
    std::vector<T> getPositions_x() const
    {
        std::vector<T> x(Positions.size(), (T)0);
        for (size_t i = 0; i < x.size(); ++i)
            x[i] = Positions[i].x;
        return x;
    }
    std::vector<T> getPositions_y() const
    {
        std::vector<T> y(Positions.size(), (T)0);
        for (size_t i = 0; i < y.size(); ++i)
            y[i] = Positions[i].y;
        return y;
    }
    std::vector<T> getVelocities_x() const
    {
        std::vector<T> Vx(Velocities.size(), (T)0);
        for (size_t i = 0; i < Vx.size(); ++i)
            Vx[i] = Velocities[i].x;
        return Vx;
    }
    std::vector<T> getVelocities_y() const
    {
        std::vector<T> Vy(Velocities.size(), (T)0);
        for (size_t i = 0; i < Vy.size(); ++i)
            Vy[i] = Velocities[i].x;
        return Vy;
    }

    // Derived types
    std::vector<T> getDisplacements_x() const
    {
        std::vector<Point<T>> Positions = getLiftedPositions();
        std::vector<T> DX(Positions.size(), (T)0);
        for (size_t i = 0; i < DX.size(); ++i)
            DX[i] = Positions[i].x - Positions[0].x;
        return DX;
    }
    std::vector<std::vector<T>> getDisplacements_x(const T &multiplier, const uint &nMoments) const
    {
        std::vector<T> DX = getDisplacements_x();
        std::vector<std::vector<T>> DXn(nMoments, std::vector<T>(DX.size(), -1));
        for (size_t i = 0; i < DX.size(); ++i)
        {
            T dxBase1 = pow(fabs(DX[i]), multiplier);
            T dxBase2 = dxBase1;
            for (uint j = 0; j < nMoments; ++j)
            {
                DXn[j][i] = dxBase1;
                dxBase1 *= dxBase2;
            }
        }
        return DXn;
    }
    std::vector<T> getTimes() const { return CumSum(getFreePaths()); }

    // Re-construct a particle object
    Particle<T> MakeParticle(const size_t &idx) const { return Particle<T>(Positions[idx], CellIDs[idx], Velocities[idx], FreePaths[idx], getTimes()[idx], ObstacleIDs[idx], SegmentIDs[idx]); }
    std::vector<Particle<T>> getParticles() const
    {
        std::vector<Particle<T>> Particles(nTotalElements, Particle<T>());
        for (size_t i = 0; i < Particles.size(); ++i)
            Particles[i] = MakeParticle(i);
        return Particles;
    }

    // Converts distretely determined particle properties to continuous ones (cannot be reversed!s)
    Trajectory<T> MakeContinuous(const std::vector<T> &Time);

    T getMeanFreePath() const;

    // Conditional trajectories
    Trajectory<T> getConditional_Space(const T &L);

    void Write(const std::string &Destination);

private:
    size_t nCurrentElements; // Current number of elements
    size_t nTotalElements;   // Expected total number at end of simulation

    // Fundamental particle types
    std::vector<Point<T>> Positions;
    std::vector<Point<cid_int>> CellIDs;
    // std::vector<Point<T>> BirkhoffPosition;
    std::vector<Point<T>> Velocities;
    std::vector<T> FreePaths;
    std::vector<int> SegmentIDs;
    std::vector<int> ObstacleIDs;
};

template <typename T>
T Trajectory<T>::getMeanFreePath() const
{
    OnlineStatistics<T> OS;
    for (size_t i = 0; i < FreePaths.size(); ++i)
        OS.Update(FreePaths[i]);
    return OS.getMean();
}

template <typename T>
Trajectory<T> Trajectory<T>::MakeContinuous(const std::vector<T> &Time)
{
    std::vector<Particle<T>> Particles = getParticles();
    Trajectory<T> ContinuousTrajectory(Time.size(), Particles[0]);
    size_t timeCount = 1;

    for (size_t i = 1; i < nTotalElements; ++i)
    {
        Particle<T> tP = Particles[i - 1]; // Particle<T> in continuous time
        while (Particles[i].TotalTime > Time[timeCount] && timeCount < Time.size())
        {
            T tau = Time[timeCount] - tP.TotalTime;
            tP.UpdateState(tau, -1, -1); // Continuous trajectory assumed to never lie on an obstacle
            ContinuousTrajectory.Push(tP);
            timeCount++;
        }
    }
    return ContinuousTrajectory;
}

// @L The location of absorbing boundary - Assumed for now to be in horiztal axis
template <typename T>
Trajectory<T> Trajectory<T>::getConditional_Space(const T &L)
{
    Trajectory<T> ConditionalTraj = this;
    std::vector<Point<T>> LiftedPositions = getLiftedPositions();
    for (size_t i = 0; i < LiftedPositions.size(); ++i)
        if (fabs(LiftedPositions[i].x) > L)
            ConditionalTraj.Pop(i + 1);
    return ConditionalTraj;
}

template <typename T>
void Trajectory<T>::Write(const std::string &Destination)
{
    std::ofstream Writer(Destination, std::ios::out);
    Writer << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;

    std::vector<Point<T>> LiftedPositions = getLiftedPositions();

    if (Writer.is_open())
    {
        size_t sz = getPositions().size();
        for (size_t i = 0; i < sz; ++i)
            Writer << LiftedPositions[i].x << " " << LiftedPositions[i].y << " " << Velocities[i].x << " " << Velocities[i].y << " " << FreePaths[i] << std::endl;
    }
    else
    {
        std::cout << "Issue opening file: " << Destination << std::endl;
    }
    Writer.close();
}

#endif