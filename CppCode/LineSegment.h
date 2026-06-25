#ifndef LINESEGMENT_H_INCLUDED
#define LINESEGMENT_H_INCLUDED

#include <iostream>
#include <iomanip>
#include "Point.h"
#include "Geometry.h"

template <typename T>
struct LineSegment
{
    LineSegment() = default;
    LineSegment(const Point<T> &p, const Point<T> &q) // : P(p), Q(q)
    {
        P = p;
        Q = q;
        InwardNormal = Point<T>(-(Q.y - P.y), (Q.x - P.x));
        InwardNormal = InwardNormal / Norm(InwardNormal);
        MidPoint = (P + Q) / static_cast<T>(2);
        Length = Distance2D<T>(P, Q);
        HalfLength = Distance2D<T>(P, MidPoint);
    };
    Point<T> P, Q;
    Point<T> InwardNormal, MidPoint;
    T Length, HalfLength;
};

template <typename T>
std::ostream &operator<<(std::ostream &os, const LineSegment<T> &L)
{
    os << std::fixed << std::setprecision(5) << L.P << " ---> " << L.Q << std::endl;
    os << "Inward normal = " << L.InwardNormal << std::endl;
    os << "Midpoint      = " << L.MidPoint << std::endl;
    os << "Length        = " << L.Length << std::endl;
    return os;
}

#endif // LINESEGMENT_H_INCLUDED
