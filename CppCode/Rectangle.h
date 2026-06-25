#ifndef RECTANGLE_H_INCLUDED
#define RECTANGLE_H_INCLUDED

#include "Polygon.h"

using fPoint = Point<FLOAT>;
class Rectangle : public Polygon
{
public:
    Rectangle(const FLOAT &width, const FLOAT &height) : Polygon({0., width, width, 0.}, {0., 0., height, height})
    {
        Width = width;
        Height = height;
    }

    Rectangle(const fPoint &start, const FLOAT &width, const FLOAT &height) : Polygon({start,
                                                                                       start + fPoint(width, 0.),
                                                                                       start + fPoint(width, height),
                                                                                       start + fPoint(0., height)})
    {
        Width = width;
        Height = height;
    }

private:
    FLOAT Width, Height;
};

#endif // RECTANGLE_H_INCLUDED
