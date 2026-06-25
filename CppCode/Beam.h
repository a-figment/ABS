#ifndef BEAM_H_INCLUDED
#define BEAM_H_INCLUDED

#include <iostream>
#include <iomanip>
#include "Point.h"
#include "Geometry.h"

// Gives the end points of a beam
template <typename T>
struct Beam
{
    Beam() = default;

    // @RayAngle The angle which relates A to B
    // @VelAngle An angle chose uniformly at random which is to be use for all particle velocities
    Beam(const Point<T> &a, const Point<T> &b, const T &rayangle, const T& velang) : A(a), B(b), RayAngle(rayangle), VelAngle(velang) //, nPointsOnBeam(npob)
    {
        Length = Distance2D<T>(A, B);
    };

    Point<T> A, B;
    T RayAngle; // The angle which relates A to B
    T VelAngle; // Velocity of particles on beam 
    T Length;
    std::vector<Point<T>> PointsOnBeam;

    // 
    Point<T> getPointOnBeam(const T& t) const {
        T x = A.x * ((T)1 - t) + B.x * t;
        T y = A.y * ((T)1 - t) + B.y * t;
        return Point<T>(x, y);
    }

    void Write(const std::string &Destination);
};


template <typename T>
void Beam<T>::Write(const std::string &Destination)
{
    std::ofstream Writer(Destination, std::ios::out | std::ios::app);
    Writer << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint << std::fixed;
    if (Writer.is_open())
        Writer << A.x << " " << A.y << " " << B.x << " " << B.y << " " << VelAngle << " " << RayAngle << " " << Length << "\n";
    else
        std::cout << "Issue opening file: " << Destination << std::endl;
    Writer.close();
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const Beam<T> &B)
{
    os << std::fixed << std::setprecision(5) << "Beam from " << B.A << " to " << B.B << " oriented with angle " << B.RayAngle * 180/3.14 << ", and velocity angle " << B.VelAngle * 180/3.14 << std::endl;
    return os;
}


#endif // BEAM_H_INCLUDED
