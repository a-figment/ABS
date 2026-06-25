#ifndef PARTICLE_H_INCLUDED
#define PARTICLE_H_INCLUDED

#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include "Geometry.h"

template <typename T>
struct Particle
{
    Particle() = default;
    Particle(const Point<T> &pos, const Point<T> &vel, const int &sid, const T &maxtime)
    {
        Pos = pos;
        Vel = vel;
        BirkhoffPos = Point<T>(-1.0, -1.0);
        CellID = Point<cid_int>(0, 0);
        FreePathLength = static_cast<T>(0);
        TotalTime = static_cast<T>(0);
        MaxTime = maxtime;
        ObstacleID = 0;
        SegmentID = sid;
        nUpdates = 0;
        onPeriodic = false;
        collisionFound = false;
    }

    Particle(const Point<T> &pos, const Point<T> &vel, const int &oid, const int &sid, const T &maxtime)
    {
        Pos = pos;
        Vel = vel;
        BirkhoffPos = Point<T>(-1.0, -1.0);
        CellID = Point<cid_int>(0, 0);
        FreePathLength = static_cast<T>(0);
        TotalTime = static_cast<T>(0);
        MaxTime = maxtime;
        ObstacleID = oid;
        SegmentID = sid;
        nUpdates = 0;
        onPeriodic = false;
        collisionFound = false;
    }

    Particle(const Point<T> &pos, const Point<T> &vel, const int &oid, const int &sid)
    {
        Pos = pos;
        Vel = vel;
        BirkhoffPos = Point<T>(-1.0, -1.0);
        CellID = Point<cid_int>(0, 0);
        FreePathLength = static_cast<T>(0);
        TotalTime = static_cast<T>(0);
        ObstacleID = oid;
        SegmentID = sid;
        nUpdates = 0;
        onPeriodic = false;
        collisionFound = false;
    }

    Particle(const Point<T> &pos, const Point<cid_int> CID, const Point<T> &vel, const T &fp, const T &tt, const int &oid, const int &sid)
    {
        Pos = pos;
        Vel = vel;
        BirkhoffPos = Point<T>(-1.0, -1.0);
        CellID = CID;
        FreePathLength = fp;
        TotalTime = tt;
        ObstacleID = oid;
        SegmentID = sid;
        nUpdates = 0;
        onPeriodic = false;
        collisionFound = false;
    }

    void UpdateState(const T &TimeToCollision, const int &OID, const int &SID)
    {
        Pos = Pos + TimeToCollision * Vel;
        BirkhoffPos = Point<T>(-1.0, -1.0);
        // CellID = Point<cid_int>(0, 0); Updated separately
        FreePathLength = TimeToCollision;
        TotalTime += TimeToCollision;
        ObstacleID = OID;
        SegmentID = SID;
        nUpdates++;
        //nUpdates += (SID != 2 || SID != 5) ? 0 : 1; // Accounts for 

        collisionFound = true;
    }

    Point<T> LPos(const Point<T>& Diameter) {
	return (Pos + Point<T>((T)CellID.x * Diameter.x, (T)CellID.y * Diameter.y)); 
    }
    void WriteState(const std::string &Destination);

    //Point<T> Diameter;            // Diameter of the billiard
    Point<T> Pos;
    Point<T> Vel;
    Point<T> BirkhoffPos;  // (s, phi) = (arc length, incidence angle)
    Point<cid_int> CellID; // Cell location of particle in plane -
    T FreePathLength;      // Length of line segment between succesive collisions
    T TotalTime;           // Sum of all free paths up to current collision
    T MaxTime;             //
    int ObstacleID;        // Many obstaces, each with 1 or more segments
    int SegmentID;         // e.g. segment of a polygon or other piecewise conncted object
    int nUpdates;          // Number of times state has been updated
    bool onPeriodic;       // Periodic boundaries
    bool collisionFound;   // Report error if false
};

template <typename T>
void Particle<T>::WriteState(const std::string &Destination)
{
    std::ofstream Writer(Destination, std::ios::out | std::ios::app);
    Writer << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint << std::fixed;
    if (Writer.is_open())
        Writer << Pos.x + CellID.x << " " << Pos.y + CellID.y << " " << Vel.x << " " << Vel.y << " " << FreePathLength << " " << ObstacleID << " " << SegmentID << "\n";
    else
        std::cout << "Issue opening file: " << Destination << std::endl;
    Writer.close();
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const Particle<T> &P)
{
    os << "Position    = " << P.Pos << "\n"
       << "Velocity    = " << P.Vel << "\n"
       << "Time        = " << P.TotalTime << "\n"
       << "Num. Update = " << P.nUpdates << "\n"
       << "Cell ID     = " << P.CellID << "\n"
       << "Obs. ID     = " << P.ObstacleID << "\n"
       << "Seg. ID     = " << P.SegmentID << "\n"
       << "On periodic = " << (P.onPeriodic ? "true" : "false") << std::endl;
    return os;
}

#endif //PARTICLE_H_INCLUDED
