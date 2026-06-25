#ifndef TRAJECTORY_H_INCLUDED
#define TRAJECTORY_H_INCLUDED

#include "Particle.h"
#include "Random.h"
#include "Ensemble.h"
// #include "Billiard.h"
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
    Trajectory(const size_t &Size, const Point<T>& diam)
    {
        nCurrentElements = 0;
        nTotalElements = Size;
        Diameter = diam;
    }

    void Reserve(const size_t &EstimatedSize, const Point<T>& diam)
    {
        nTotalElements = 0;
        nCurrentElements = 0;

        LiftedPositions.reserve(EstimatedSize);
        Positions.reserve(EstimatedSize);
        CellIDs.reserve(EstimatedSize);
        Velocities.reserve(EstimatedSize);
        FreePaths.reserve(EstimatedSize);
        Time.reserve(EstimatedSize);
        SegmentIDs.reserve(EstimatedSize);
        ObstacleIDs.reserve(EstimatedSize);

        Diameter = diam;
    }

    Trajectory(const size_t &Size, const Particle<T> &IP, const Point<T>& diam)
    {
        nCurrentElements = 0;
        nTotalElements = Size;

        LiftedPositions = std::vector<Point<T>>(Size, Point<T>());
        Positions = std::vector<Point<T>>(Size, Point<T>());
        BirkhoffPositions = std::vector<Point<T>>(Size, Point<T>());
        CellIDs = std::vector<Point<cid_int>>(Size, Point<cid_int>());
        Velocities = std::vector<Point<T>>(Size, Point<T>());
        FreePaths = std::vector<T>(Size, (T)0);
        Time = std::vector<T>(Size, (T)0);
        SegmentIDs = std::vector<int>(Size, -2); // -1 => particle is floating. -2 => it is undefined
        ObstacleIDs = std::vector<int>(Size, -2);

        Diameter = diam;
        Push(IP);
    }

    void Push(const Particle<T> &P)
    {
        LiftedPositions[nCurrentElements] = P.Pos + Point<T>((T)P.CellID.x * Diameter.x, (T)P.CellID.y * Diameter.y); 
	//std::cout << "P = " << P.Pos << std::endl;
	//std::cout << "p = " << Point<T>((T)P.CellID.x * Diameter.x, (T)P.CellID.y * Diameter.y) << std::endl;
        Positions[nCurrentElements] = P.Pos;
        BirkhoffPositions[nCurrentElements] = P.BirkhoffPos;
        CellIDs[nCurrentElements] = P.CellID;
        Velocities[nCurrentElements] = P.Vel;
        FreePaths[nCurrentElements] = P.FreePathLength;
        Time[nCurrentElements] = P.TotalTime;
        SegmentIDs[nCurrentElements] = P.SegmentID;
        ObstacleIDs[nCurrentElements] = P.ObstacleID;
        nCurrentElements++;
    }

    void Push_(const Particle<T> &P)
    {
        LiftedPositions.push_back(P.Pos + Point<T>((T)P.CellID.x * Diameter.x, (T)P.CellID.y * Diameter.y)); // Lifted position
        Positions.push_back(P.Pos);
        BirkhoffPositions.push_back(P.BirkhoffPos);
        CellIDs.push_back(P.CellID);
        Velocities.push_back(P.Vel);
        FreePaths.push_back(P.FreePathLength);
        Time.push_back(P.TotalTime);
        SegmentIDs.push_back(P.SegmentID);
        ObstacleIDs.push_back(P.ObstacleID);
        nCurrentElements++;
        // nTotalElements = (nCurrentElements > nTotalElements ? nCurrentElements : nTotalElements);
    }

    void Pop(const size_t &offset)
    {
        LiftedPositions.erase(Positions.begin() + offset, Positions.end());
        Positions.erase(Positions.begin() + offset, Positions.end());
        BirkhoffPositions.erase(BirkhoffPositions.begin() + offset, BirkhoffPositions.end());
        CellIDs.erase(CellIDs.begin() + offset, CellIDs.end());
        Velocities.erase(Velocities.begin() + offset, Velocities.end());
        FreePaths.erase(FreePaths.begin() + offset, FreePaths.end());
        Time.erase(Time.begin() + offset, Time.end());
        SegmentIDs.erase(SegmentIDs.begin() + offset, SegmentIDs.end());
        ObstacleIDs.erase(ObstacleIDs.begin() + offset, ObstacleIDs.end());
        nTotalElements = Positions.size();
    }

    // Base types
    size_t getNumElements() const { return nCurrentElements; }
    size_t getTotalLength() const { return (nCurrentElements > nTotalElements ? nCurrentElements : nTotalElements); }

    std::vector<Point<T>> getPositions() const { return Positions; }
    std::vector<Point<T>> getBirkhoffPositions() const { return BirkhoffPositions; }

    std::vector<Point<T>> getLiftedPositions() const
    {
        return LiftedPositions;
    }
    std::vector<Point<cid_int>> getCellIDs() const { return CellIDs; }

    std::vector<cid_int> getCellIDs_x() const { 
	std::vector<cid_int> cidx(CellIDs.size(), 0);
	for (size_t i = 0; i < CellIDs.size(); ++i)
	    cidx[i] = CellIDs[i].x;
	return cidx; 
    }
    std::vector<Point<T>> getVelocities() const { return Velocities; }
    std::vector<T> getFreePaths() const { return FreePaths; }
    std::vector<T> getTime() const { return Time; }
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

    std::vector<T> getArcLengths() const
    {
        std::vector<T> s(BirkhoffPositions.size(), (T)0);
        for (size_t i = 0; i < s.size(); ++i)
            s[i] = BirkhoffPositions[i].x;
        return s;
    }

    std::vector<T> getIncidenceAngles() const
    {
        std::vector<T> phi(BirkhoffPositions.size(), (T)0);
        for (size_t i = 0; i < phi.size(); ++i)
            phi[i] = BirkhoffPositions[i].y;
        return phi;
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
            Vy[i] = Velocities[i].y;
        return Vy;
    }
    std::vector<T> getVelocities_a() const
    {
        std::vector<T> Vx = getVelocities_x(); 
        std::vector<T> Vy = getVelocities_y(); 
        std::vector<T> Angle(Velocities.size()); 
        for (size_t i = 0; i < Vy.size(); ++i)
            Angle[i] = atan2(Velocities[i].y, Velocities[i].x);
        return Angle;
    }


    // Distance between each collision on lifted domain
    std::vector<T> getLiftedFreePaths() const
    {
        std::vector<Point<T>> Positions = getLiftedPositions();
        std::vector<T> R(Positions.size(), (T)0); // n collisions, n-1 free paths, Position[0] occurs at t = 0
        for (size_t i = 1; i < R.size(); ++i)
            R[i] = Distance2D<T>(Positions[i-1], Positions[i]);
        return R;
    }


    // Derived types
    std::vector<T> getLiftedPositions_x() const
    {
        std::vector<Point<T>> Positions = getLiftedPositions();
        std::vector<T> X(Positions.size(), (T)0);
        for (size_t i = 0; i < X.size(); ++i)
            X[i] = Positions[i].x;
        return X;
    }

    std::vector<T> getLiftedPositions_y() const
    {
        std::vector<Point<T>> Positions = getLiftedPositions();
        std::vector<T> Y(Positions.size(), (T)0);
        for (size_t i = 0; i < Y.size(); ++i)
            Y[i] = Positions[i].y;
        return Y;
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

    std::vector<T> getDisplacements_y() const
    {
        std::vector<Point<T>> Positions = getLiftedPositions();
        std::vector<T> DY(Positions.size(), (T)0);
        for (size_t i = 0; i < DY.size(); ++i)
            DY[i] = Positions[i].y - Positions[0].y;
        return DY;
    }

    std::vector<T> getDisplacements_xy() const
    {
        std::vector<Point<T>> Positions = getLiftedPositions();
        std::vector<T> R(Positions.size(), (T)0);
        for (size_t i = 0; i < R.size(); ++i)
            R[i] = (Positions[i].x - Positions[0].x);// + (Positions[i].y - Positions[0].y);
        return R;
    }

    std::vector<T> getVAC_x() const
    {
        std::vector<T> Vx = getVelocities_x();
        std::vector<T> Corr(Vx.size(), (T)0);
        for (size_t i = 0; i < Corr.size(); ++i)
            Corr[i] = Vx[i] * Vx[0];
        return Corr;
    }

    std::vector<T> getBirkVAC_x() const
    {
        std::vector<T> Inc = getIncidenceAngles();
        std::vector<T> Corr(Inc.size(), (T)0);
        for (size_t i = 1; i < Corr.size(); ++i)
            Corr[i] = cos(Inc[i]) * cos(Inc[1]);
        return Corr;
    }

    std::vector<T> getVAC_xy() const
    {
        std::vector<Point<T>> Vxy = getVelocities();
        Point<T> Vxy0 = Vxy[0];
        std::vector<T> Corr(Vxy.size(), (T)0);
        for (size_t i = 0; i < Corr.size(); ++i)
            Corr[i] = Dot<T>(Vxy0, Vxy[i]);
        return Corr;
    }

    std::vector<T> getIVAC_x() const
    {
        T Vx0 = getVelocities_x()[0];
        std::vector<T> DX = getDisplacements_x();
        std::vector<T> ICorr(DX.size(), (T)0);

        for (size_t i = 0; i < ICorr.size(); ++i)
            ICorr[i] = DX[i] * Vx0;
        return ICorr;
    }

    std::vector<T> getDisplacements_abs_x() const
    {
        std::vector<Point<T>> Positions = getLiftedPositions();
        std::vector<T> DX(Positions.size(), (T)0);
        for (size_t i = 0; i < DX.size(); ++i)
            DX[i] = fabs(Positions[i].x - Positions[0].x);
        return DX;
    }

    std::vector<std::vector<T>> getDisplacements_x(const T &multiplier, const uint &nMoments) const
    {
        std::vector<T> DX = getDisplacements_x();
        std::vector<std::vector<T>> DXn(nMoments, std::vector<T>(DX.size(), -1));
        for (size_t i = 0; i < DX.size(); ++i)
        {
            //T dxBase1 = pow(std::abs(DX[i]), multiplier);
            T dxBase1 = pow(DX[i], multiplier);
            T dxBase2 = dxBase1;
            for (size_t j = 0; j < nMoments; ++j)
            {
                DXn[j][i] = dxBase1;
                dxBase1 *= dxBase2;
            }
        }
        return DXn;
    }

    std::vector<std::vector<T>> getDisplacementPowers_x(const std::vector<T> &moments) const
    {
	std::vector<T> DX = getDisplacements_x();
	std::vector<std::vector<T>> DXn(moments.size(), std::vector<T>(DX.size(), -1));

	for (size_t i = 0; i < moments.size(); ++i) {
	    for (size_t j = 0; j < DX.size(); ++j) {
		DXn[i][j] = pow(std::abs(DX[j]), moments[i]);
	    }
	}
	return DXn;
    }



    std::vector<std::vector<T>> getDisplacements_xy(const T &multiplier, const uint &nMoments) const
    {
        std::vector<T> DX = getDisplacements_x();
        std::vector<T> DY = getDisplacements_y();

        std::vector<std::vector<T>> DXn(nMoments, std::vector<T>(DX.size(), -1));
        for (size_t i = 0; i < DX.size(); ++i)
        {
            T R = sqrt(DX[i] * DX[i] + DY[i] * DY[i]);
            T dxBase1 = pow(R, multiplier);
            T dxBase2 = dxBase1;
            for (uint j = 0; j < nMoments; ++j)
            {
                DXn[j][i] = dxBase1;
                dxBase1 *= dxBase2;
            }
        }
        return DXn;
    }

    // Re-construct a particle object
    Particle<T> MakeParticle(const size_t &idx) const { return Particle<T>(Positions[idx], CellIDs[idx], Velocities[idx], FreePaths[idx], Time[idx], ObstacleIDs[idx], SegmentIDs[idx]); }

    std::vector<Particle<T>> getParticles() const
    {
        std::vector<Particle<T>> Particles(Time.size());
        for (size_t i = 0; i < Particles.size(); ++i)
            Particles[i] = MakeParticle(i);
        return Particles;
    }

    // Converts distretely determined particle properties to continuous ones (cannot be reversed!s)
    Trajectory<T> MakeContinuous(const std::vector<T> &Time);

    T getMeanFreePath() const;

    // Conditional trajectories
    Trajectory<T> getConditional_Space(const T &L);

    void WritePositions_App(const std::string &Destination);
    void Write(const std::string &Destination);
    void WriteVelocities(const std::string &Destination);
    void WritePositions(const std::string &Destination);
    size_t nCollisions;
private:
    size_t nCurrentElements; // Current number of elements
    size_t nTotalElements;   // Expected total number at end of simulation

    // Fundamental particle types
    std::vector<Point<T>> Positions;
    std::vector<Point<cid_int>> CellIDs;
    std::vector<Point<T>> BirkhoffPositions;
    std::vector<Point<T>> Velocities;
    std::vector<T> FreePaths;
    std::vector<T> Time;
    std::vector<int> SegmentIDs;
    std::vector<int> ObstacleIDs;

    std::vector<Point<T>> LiftedPositions;
    Point<T> Diameter;

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
    // Trajectory<T> ContinuousTrajectory(Time.size(), Particles[0]);
    Trajectory<T> ContinuousTrajectory;
    ContinuousTrajectory.Reserve(Time.size());
    ContinuousTrajectory.Push_(Particles[0]);

    size_t timeCount = 1;

    for (size_t i = 1; i < Particles.size(); ++i)
    {
        Particle<T> tP = Particles[i - 1]; // Particle<T> in continuous time
        while (Particles[i].TotalTime > Time[timeCount] && timeCount < Time.size())
        {
            T tau = Time[timeCount] - tP.TotalTime;
            tP.UpdateState(tau, -1, -1); // Continuous trajectory assumed to never lie on an obstacle
            ContinuousTrajectory.Push_(tP);
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
            Writer << LiftedPositions[i].x << " " << LiftedPositions[i].y << " " << Velocities[i].x << " " << Velocities[i].y << " "
                   << BirkhoffPositions[i].x << " " << BirkhoffPositions[i].y << " " << FreePaths[i] << std::endl;
    }
    else
    {
        std::cout << "Issue opening file: " << Destination << std::endl;
    }
    Writer.close();
}

template <typename T>
void Trajectory<T>::WritePositions(const std::string &Destination)
{
    std::ofstream Writer(Destination, std::ios::out);
    Writer << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;

    std::vector<Point<T>> LiftedPositions = getLiftedPositions();

    if (Writer.is_open())
    {
        size_t sz = getPositions().size();
        for (size_t i = 0; i < sz; ++i)
            Writer << LiftedPositions[i].x << " " << LiftedPositions[i].y << std::endl;
    }
    else
    {
        std::cout << "Issue opening file: " << Destination << std::endl;
    }
    Writer.close();
}

// Only to be used for debugging - very inefficient. Make a Writer class already!!
template <typename T>
void Trajectory<T>::WritePositions_App(const std::string &Destination)
{
    std::ofstream Writer(Destination, std::ios::out | std::ios::app);
    Writer << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;

    std::vector<Point<T>> LiftedPositions = getLiftedPositions();

    if (Writer.is_open())
    {
        size_t sz = getPositions().size();
        for (size_t i = 0; i < sz; ++i)
            Writer << LiftedPositions[i].x << " " << LiftedPositions[i].y << std::endl;
    }
    else
    {
        std::cout << "Issue opening file: " << Destination << std::endl;
    }
    Writer.close();
}


template <typename T>
void Trajectory<T>::WriteVelocities(const std::string &Destination)
{
    std::ofstream Writer(Destination, std::ios::out);
    Writer << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;

    std::vector<Point<T>> LiftedPositions = getLiftedPositions();

    if (Writer.is_open())
    {
        size_t sz = getPositions().size();
        for (size_t i = 0; i < sz; ++i)
            Writer << Velocities[i].x << " " << Velocities[i].y << std::endl;
    }
    else
    {
        std::cout << "Issue opening file: " << Destination << std::endl;
    }
    Writer.close();
}


#endif

// #ifndef TRAJECTORY_H_INCLUDED
// #define TRAJECTORY_H_INCLUDED

// #include "Particle.h"
// #include "Random.h"
// #include "Ensemble.h"
// #include <sstream>

// template <typename T>
// struct TrajectoryStatistics
// {

//     TrajectoryStatistics()
//     {
//         MeanNumberOfCollisions = 0;
//         MeanTimeAferCollisions = 0;
//         MFP = 0;
//     }
//     T MeanNumberOfCollisions;
//     T MeanTimeAferCollisions;
//     T MFP;
// };

// template <typename T>
// struct EscapeData
// {
//     EscapeData(const std::vector<T> &nEscapeEvents)
//     {
//         EscapeTimes.resize(nEscapeEvents.size());
//         EscapePosition_x.resize(nEscapeEvents.size());
//         EscapePosition_y.resize(nEscapeEvents.size());
//         EscapeDirection.resize(nEscapeEvents.size());
//         xPositionSum.resize(nEscapeEvents.size());
//         nEscapers.resize(nEscapeEvents.size());
//     }

//     EscapeData(const size_t &sz)
//     {
//         EscapeTimes.resize(sz);
//         EscapePosition_x.resize(sz);
//         EscapePosition_y.resize(sz);
//         EscapeDirection.resize(sz);
//         xPositionSum.resize(sz);
//         nEscapers.resize(sz);
//     }

//     std::vector<T> EscapeTimes;
//     std::vector<T> EscapePosition_x;
//     std::vector<T> EscapePosition_y;
//     std::vector<T> EscapeDirection;
//     std::vector<T> xPositionSum; // Sum of x-positions up until an escape
//     std::vector<int> nEscapers;  // Sum of x-positions up until an escape
// };

// template <typename T>
// class Trajectory
// {
// public:
//     Trajectory()
//     {
//         nElems = 0;
//         Length = 0;
//     }

//     // @Length - Total number of particles to consider
//     Trajectory(const size_t &length)
//     {
//         nElems = 0;
//         Length = length;
//         Particles.resize(Length);
//         for (size_t i = 0; i < length; ++i)
//             Particles[i] = Particle<T>();
//     }

//     Trajectory(const size_t &length, const Particle<T> &IP)
//     {
//         nElems = 0;
//         Length = length;
//         Particles.resize(Length);
//         for (size_t i = 0; i < length; ++i)
//             Particles[i] = Particle<T>();
//         Push(IP);
//     }

//     void Reserve(const size_t &EstimatedSize) { Particles.reserve(EstimatedSize); }
//     void Push(const Particle<T> &P) { Particles[nElems++] = P; }
//     void Pop(const size_t &offset)
//     {
//         Particles.erase(Particles.begin() + offset, Particles.end());
//         Length = Particles.size();
//     }
//     void Push_(const Particle<T> &P)
//     {
//         Particles.push_back(P);
//         nElems++;
//     }

//     std::vector<Point<T>> getPositions() const
//     {
//         std::vector<Point<T>> Positions(Particles.size(), Point<T>());
//         for (size_t i = 0; i < Particles.size(); ++i)
//             Positions[i] = LiftPosition(Particles[i]);
//         return Positions;
//     }
//     std::vector<T> getDisplacements_x() const
//     {
//         std::vector<Point<T>> LiftedPositions = getPositions();
//         std::vector<T> DX(Particles.size());
//         for (size_t i = 0; i < Particles.size(); ++i)
//             DX[i] = LiftedPositions[i].x - LiftedPositions[0].x;
//         return DX;
//     }
//     std::vector<T> getDisplacementsUpToL(const T &L) const
//     {
//         std::vector<Point<T>> LiftedPositions = getPositions();
//         std::vector<T> DX(Particles.size());
//         for (size_t i = 0; i < Particles.size(); ++i)
//         {
//             DX[i] = LiftedPositions[i].x - LiftedPositions[0].x;
//             if (fabs(LiftedPositions[i].x) > L)
//                 break;
//         }
//         return DX;
//     }
//     std::vector<T> getDisplacements_x(const T &power) const
//     {
//         std::vector<Point<T>> LiftedPositions = getPositions();
//         std::vector<T> DXn(Particles.size());
//         for (size_t i = 0; i < Particles.size(); ++i)
//             DXn[i] = pow(LiftedPositions[i].x - LiftedPositions[0].x, power);
//         return DXn;
//     }
//     std::vector<std::vector<T>> getDisplacements_x(const T &multiplier, const uint &nMoments) const
//     {
//         std::vector<T> DX = getDisplacements_x();
//         std::vector<std::vector<T>> DXn(nMoments, std::vector<T>(Particles.size(), -1));

//         for (size_t i = 0; i < DX.size(); ++i)
//         {
//             T dxBase1 = pow(fabs(DX[i]), multiplier);
//             T dxBase2 = dxBase1;

//             for (uint j = 0; j < nMoments; ++j)
//             {
//                 DXn[j][i] = dxBase1;
//                 dxBase1 *= dxBase2;
//             }
//         }

//         return DXn;
//     }
//     std::vector<T> getDisplacements_x2() const
//     {
//         std::vector<Point<T>> LiftedPositions = getPositions();
//         std::vector<T> DX2(Particles.size());
//         for (size_t i = 0; i < Particles.size(); ++i)
//             DX2[i] = (LiftedPositions[i].x - LiftedPositions[0].x) * (LiftedPositions[i].x - LiftedPositions[0].x);
//         return DX2;
//     }
//     std::vector<T> getTimes() const
//     {
//         std::vector<T> Time(Particles.size());
//         for (size_t i = 0; i < Particles.size(); ++i)
//             Time[i] = Particles[i].TotalTime;
//         return Time;
//     }
//     std::vector<T> getFreePaths() const
//     {
//         std::vector<T> FreePaths(Particles.size());
//         for (size_t i = 0; i < Particles.size(); ++i)
//             FreePaths[i] = Particles[i].FreePathLength;
//         return FreePaths;
//     }
//     std::vector<int> getSegmentIDs() const
//     {
//         std::vector<int> IDs(Particles.size());
//         for (size_t i = 0; i < Particles.size(); ++i)
//             IDs[i] = Particles[i].SegmentID;
//         return IDs;
//     }
//     std::vector<T> getVelocities_x() const
//     {
//         std::vector<T> xVel(Particles.size());
//         for (size_t i = 0; i < Particles.size(); ++i)
//             xVel[i] = Particles[i].Vel.x;
//         return xVel;
//     }

//     size_t getNumElements() const { return nElems; }
//     size_t getLength() const { return Particles.size(); }
//     std::vector<Particle<T>> getParticles() const { return Particles; }

//     void MakeContinuous(const std::vector<T> &Time);

//     Particle<T> getParticle(const size_t &index) const { return Particles[index]; }
//     std::vector<Point<T>> getEscapePositions(const std::vector<T> &L);
//     std::vector<T> getEscapeTimes(const std::vector<T> &L);
//     std::vector<T> getEscapeTimes(const std::vector<std::string> &Symb);
//     std::vector<T> getEscapeTimesAtRandomCollision(const int &nToGet);
//     std::vector<T> getEscapeDirection(const std::vector<T> &L);

//     std::vector<T> getPositionsUpToEscape(const T &L);
//     std::vector<T> getDiffsUpToEscape(const T &L);

//     EscapeData<T> getEscapeData(const std::vector<T> &L);
//     EscapeData<T> getEscapeData(const std::vector<std::string> &Symbols);

//     T getMeanFreePath() const;

//     void Write(const std::string &Destination);

// private:
//     Point<T> LiftPosition(const Particle<T> &P) const;
//     size_t nElems;                      // Current number of elements
//     size_t Length;                      // Trajectory length (total number of elements)
//     std::vector<Particle<T>> Particles; // Forms a trajectory
// };

// template <typename T>
// Point<T> Trajectory<T>::LiftPosition(const Particle<T> &P) const
// {
//     T xLift = P.Pos.x + static_cast<T>(P.CellID.x);
//     T yLift = P.Pos.y + static_cast<T>(P.CellID.y);
//     return Point<T>(xLift, yLift);
// }

// template <typename T>
// EscapeData<T> Trajectory<T>::getEscapeData(const std::vector<T> &L)
// {
//     std::vector<Point<T>> Positions = getPositions();
//     EscapeData<T> ED(L);
//     size_t startidx = 0;

//     for (size_t i = 0; i < L.size(); ++i)
//     {
//         for (size_t j = startidx; j < Positions.size(); ++j)
//         {
//             if (fabs(Positions[j].x) > L[i]) // Escape occurred
//             {
//                 Particle<T> P = Particles[j - 1]; // Particle prior to esacpe
//                 T l = L[i];
//                 if (P.Vel.x < static_cast<T>(0))
//                     l = -L[i];
//                 Point<T> EPos = getPointOnVerticalLine(Positions[j - 1], P.Vel.y / P.Vel.x, l);
//                 const T d = Distance2D<T>(P.Pos, EPos);

//                 // Escape data
//                 ED.EscapeTimes[i] = P.TotalTime + d;
//                 ED.EscapePosition_x[i] = EPos.x;
//                 ED.EscapePosition_y[i] = EPos.y;
//                 ED.EscapeDirection[i] = atan2(P.Vel.y, P.Vel.x);

//                 // std::vector<T>::iterator upTo{Positions.begin() + j - 1};

//                 // ED.xPositionSum[i] = std::accumulate(Positions.begin(), Positions.begin() + (j - 1), (T)0) + L[i];
//                 // for (auto x : Positions)
//                 // ED.xPositionSum[i] += abs(x.x);
//                 // std::transform_reduce(Positions.begin(), vec.end(),0,std::plus<>{}, static_cast<double (*)(double)>(std::fabs));
//                 ED.nEscapers[i]++;

//                 startidx = j - 1;
//                 break;
//             }
//         }
//     }
//     return ED;
// }

// template <typename T>
// std::vector<T> Trajectory<T>::getPositionsUpToEscape(const T &L)
// {
//     std::vector<Point<T>> Positions = getPositions();
//     size_t startidx = 0;

//     std::vector<T> xPositions;
//     xPositions.reserve(Positions.size());

//     for (size_t j = startidx; j < Positions.size(); ++j)
//     {
//         if (fabs(Positions[j].x) > L) // Escape occurred
//         {
//             for (size_t k = 0; k < j; ++k)
//                 xPositions.push_back(Positions[k].x);
//             startidx = j - 1;
//             break;
//         }
//     }
//     return xPositions;
// }

// template <typename T>
// std::vector<T> Trajectory<T>::getDiffsUpToEscape(const T &L)
// {
//     std::vector<Point<T>> Positions = getPositions();
//     size_t startidx = 0;

//     std::vector<T> xDiffs;
//     xDiffs.reserve(Positions.size());

//     for (size_t j = startidx; j < Positions.size(); ++j)
//     {
//         if (fabs(Positions[j].x) > L) // Escape occurred
//         {
//             for (size_t k = 1; k < j; ++k)
//                 xDiffs.push_back(Positions[k].x - Positions[k - 1].x);
//             startidx = j - 1;
//             break;
//         }
//     }
//     return xDiffs;
// }

// template <typename T>
// EscapeData<T> Trajectory<T>::getEscapeData(const std::vector<std::string> &Symbols)
// {
//     EscapeData<T> ED(Symbols.size());
//     std::vector<int> IDs = getSegmentIDs();

//     // Convert sequence to string
//     std::ostringstream oss;
//     std::copy(IDs.begin() + 1, IDs.end(), std::ostream_iterator<int>(oss, ""));
//     std::string sID = oss.str();

//     for (size_t i = 0; i < Symbols.size(); ++i)
//     {
//         size_t idx = sID.find(Symbols[i], 0);
//         if (idx != std::string::npos)
//         {
//             ED.EscapeTimes[i] = Particles[idx].TotalTime;
//             ED.EscapePosition_x[i] = Particles[idx].Pos.x;
//             ED.EscapePosition_y[i] = Particles[idx].Pos.y;
//             ED.EscapeDirection[i] = atan2(Particles[idx].Vel.y, Particles[idx].Vel.x);
//         }
//         else
//         {
//             ED.nEscapers[i]++;
//         }
//     }
//     return ED;
// }

// // @brief This can only work is particles are collision particles
// template <typename T>
// std::vector<T> Trajectory<T>::getEscapeTimes(const std::vector<T> &L)
// {
//     std::vector<Point<T>> Positions = getPositions();
//     std::vector<T> EscapeTimes(L.size(), INFINITY);
//     size_t startidx = 0;

//     for (size_t i = 0; i < L.size(); ++i)
//     {
//         for (size_t j = startidx; j < Positions.size(); ++j)
//         {
//             if (fabs(Positions[j].x) > L[i]) // Escape occurred
//             {
//                 if (j == 0)
//                     int x = 0;
//                 Particle<T> P = Particles[j - 1]; // Particle prior to escape
//                 T l = L[i];
//                 if (P.Vel.x < static_cast<T>(0))
//                     l = -L[i];
//                 Point<T> EPos = getPointOnVerticalLine(Positions[j - 1], P.Vel.y / P.Vel.x, l);
//                 const T d = Distance2D<T>(P.Pos, EPos);
//                 EscapeTimes[i] = P.TotalTime + d;
//                 startidx = j - 1;
//                 break;
//             }
//         }
//     }
//     return EscapeTimes;
// }

// template <typename T>
// std::vector<T> Trajectory<T>::getEscapeTimes(const std::vector<std::string> &Symb)
// {
//     std::vector<Point<T>> Positions = getPositions();
//     std::vector<int> IDs = getSegmentIDs();

//     // Convert sequence to string
//     std::ostringstream oss;
//     std::copy(IDs.begin() + 1, IDs.end(), std::ostream_iterator<int>(oss, ""));
//     std::string sID = oss.str();

//     std::vector<T> EscapeTimes(Symb.size(), INFINITY);

//     for (size_t i = 0; i < Symb.size(); ++i)
//     {
//         size_t idx = sID.find(Symb[i], 0);
//         if (idx != std::string::npos)
//             EscapeTimes[i] = Particles[idx].TotalTime;
//     }
//     return EscapeTimes;
// }

// template <typename T>
// std::vector<Point<T>> Trajectory<T>::getEscapePositions(const std::vector<T> &L)
// {
//     std::vector<Point<T>> Positions = getPositions();
//     std::vector<Point<T>> EPos(L.size());
//     size_t startidx = 0;

//     for (size_t i = 0; i < L.size(); ++i)
//     {
//         for (size_t j = startidx; j < Positions.size(); ++j)
//         {
//             if (fabs(Positions[j].x) > L[i]) // Escape occurred
//             {
//                 Particle<T> P = Particles[j - 1]; // Particle prior to escape
//                 T l = L[i];
//                 if (P.Vel.x < static_cast<T>(0))
//                     l = -L[i];
//                 EPos[i] = getPointOnVerticalLine(Positions[j - 1], P.Vel.y / P.Vel.x, l);
//                 startidx = j - 1;
//                 break;
//             }
//         }
//     }
//     return EPos;
// }

// template <typename T>
// std::vector<T> Trajectory<T>::getEscapeDirection(const std::vector<T> &L)
// {
//     std::vector<Point<T>> Positions = getPositions();
//     std::vector<T> Direction(L.size(), -1);
//     size_t startidx = 0;

//     for (size_t i = 0; i < L.size(); ++i)
//     {
//         for (size_t j = startidx; j < Positions.size(); ++j)
//         {
//             if (fabs(Positions[j].x) > L[i]) // Escape occurred
//             {
//                 Particle<T> P = Particles[j - 1]; // Particle prior to escape
//                 Direction[i] = atan2(P.Vel.y, P.Vel.x);
//                 startidx = j - 1;
//                 break;
//             }
//         }
//     }
//     return Direction;
// }

// template <typename T>
// std::vector<T> Trajectory<T>::getEscapeTimesAtRandomCollision(const int &nToGet)
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
// }

// template <typename T>
// T Trajectory<T>::getMeanFreePath() const
// {
//     OnlineStatistics<T> OS;
//     for (size_t i = 0; i < Particles.size(); ++i)
//         OS.Update(Particles[i].FreePathLength);
//     return OS.getMean();
// }

// template <typename T>
// void Trajectory<T>::MakeContinuous(const std::vector<T> &Time)
// {
//     std::vector<Particle<T>> contParticles;
//     contParticles.reserve(Particles.size());
//     contParticles.push_back(Particles[0]);
//     size_t timeCount = 1;

//     for (size_t i = 1; i < Particles.size(); ++i)
//     {
//         Particle<T> tP = Particles[i - 1]; // Particle<T> in continuous time
//         while (Particles[i].TotalTime > Time[timeCount] && timeCount < Time.size())
//         {
//             T tau = Time[timeCount] - tP.TotalTime;
//             tP.UpdateState(tau, -1, -1); // Continuous trajectory assumed to never lie on an obstacle
//             contParticles.push_back(tP);
//             timeCount++;
//         }
//     }
//     Particles = contParticles;
// }

// template <typename T>
// void Trajectory<T>::Write(const std::string &Destination)
// {
//     std::ofstream Writer(Destination, std::ios::out);
//     Writer << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;
//     std::vector<Point<T>> Positions = getPositions();
//     std::vector<T> FP = getFreePaths();
//     if (Writer.is_open())
//     {
//         size_t sz = getPositions().size();
//         for (size_t i = 0; i < sz; ++i)
//             Writer << Positions[i].x << " " << Positions[i].y << " " << FP[i] << std::endl;
//     }
//     else
//     {
//         std::cout << "Issue opening file: " << Destination << std::endl;
//     }
//     Writer.close();
// }

// #endif

// // #ifndef TRAJECTORY_H_INCLUDED
// // #define TRAJECTORY_H_INCLUDED

// // #include "Particle.h"

// // template <typename T>
// // class Trajectory
// // {
// // public:
// //     Trajectory()
// //     {
// //         nElems = 0;
// //         Length = 0;
// //     }

// //     // @Length - Number of collisions, total time units, ...
// //     Trajectory(const size_t &length)
// //     {
// //         nElems = 0;
// //         Length = length;
// //         xPositions.resize(Length);
// //         yPositions.resize(Length);
// //         FreePaths.resize(Length);
// //     }

// //     void Append(const Trajectory<T> &Traj)
// //     {
// //         std::vector<Point<T>> p = Traj.getPositions();
// //         std::vector<T> tau = Traj.getFreePaths();
// //         size_t n = getNumElements();
// //         for (size_t i = 0; i < Traj.getNumElements(); ++i)
// //         {
// //             xPositions[i + n] = p[i].x;
// //             yPositions[i + n] = p[i].y;
// //             FreePaths[i + n] = tau[i];
// //             nElems++;
// //         }
// //     }

// //     void Push(const Particle<T> &P)
// //     {
// //         xPositions[nElems] = P.Pos.x + static_cast<T>(P.CellID.x);
// //         yPositions[nElems] = P.Pos.y;
// //         FreePaths[nElems] = P.FreePathLength;
// //         nElems++;
// //     }

// //     void Push_(const Particle<T> &P)
// //     {
// //         xPositions.push_back(P.Pos.x + static_cast<T>(P.CellID.x));
// //         yPositions.push_back(P.Pos.y);
// //         FreePaths.push_back(P.FreePathLength);
// //         nElems++;
// //     }

// //     std::vector<Point<T>> getPositions() const
// //     {
// //         std::vector<Point<T>> Positions(xPositions.size());
// //         for (size_t i = 0; i < xPositions.size(); ++i)
// //             Positions[i] = Point<T>(xPositions[i], yPositions[i]);
// //         return Positions;
// //     }

// //     std::vector<T> getDisplacements_x2() const
// //     {
// //         std::vector<T> DX(xPositions.size());
// //         for (size_t i = 0; i < xPositions.size(); ++i)
// //             DX[i] = (xPositions[i] - xPositions[0]) * (xPositions[i] - xPositions[0]);
// //         return DX;
// //     }

// //     std::vector<T> getTimes() const { return CumSum(FreePaths); }
// //     std::vector<T> getFreePaths() const { return FreePaths; }
// //     size_t getNumElements() const { return nElems; }

// //     void Write(const std::string &Destination)
// //     {
// //         std::ofstream Writer(Destination, std::ios::out);
// //         Writer << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;
// //         if (Writer.is_open())
// //         {
// //             for (size_t i = 0; i < xPositions.size(); ++i)
// //                 Writer << xPositions[i] << " " << yPositions[i] << " " << FreePaths[i] << "\n";
// //         }
// //         else
// //         {
// //             std::cout << "Issue opening file: " << Destination << std::endl;
// //         }
// //         Writer.close();
// //     }

// // private:
// //     size_t nElems; // Current number of elements
// //     size_t Length; // Trajectory length (indices)
// //     std::vector<T> xPositions;
// //     std::vector<T> yPositions;
// //     std::vector<T> FreePaths; // Time between two collisions
// // };

// // #endif
