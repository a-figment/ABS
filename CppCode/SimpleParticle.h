#ifndef SIMPLETRAJECTORY_H_INCLUDED
#define SIMPLETRAJECTORY_H_INCLUDED

#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include "Geometry.h"

template <typename T>
struct SimpleTrajectory
{
    Particle(const Point<T> &initial_pos, )
    {

    }


    //Point<T> Diameter;            // Diameter of the billiard
    Point<T> Pos;
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


#endif //PARTICLE_H_INCLUDED
