#ifndef CIRCLE_H_INCLUDED
#define CIRCLE_H_INCLUDED

#include <iomanip>
#include <vector>
#include "Geometry.h"
#include "Particle.h"
#include "Obstacle.h"

// TODO - The whole class!
template <typename T>
class Circle : public Obstacle<T>
{

public:
    // Circle(const T& radius, const Point<T>& center) : Radius(radius), Center(center) {}
    Circle(const T &radius, const Point<T> &center) : Obstacle<T>()
    {
        Radius = radius;
        Center = center;

        setArea();
        setPerimeter();
        setCumulativePerimeter();
        setBoundingBox();
    }

    virtual T getArcDistance(const Point<T> &P, const int &ID)
    {
        return 0;
    }

    virtual T getIncidenceAngle(const Point<T> &V, const int &ID)
    {
        return 0;
    }

    virtual void setCumulativePerimeter()
    {
        this->CumulativePerimeter = {static_cast<T>(0)}; // Arc length starts from (Radius + Center.x,0)
    }

    virtual void setPeriodicIDs(const std::vector<int> &ids)
    {
        this->PeriodicIDs[0] = std::make_pair(-1, -1); // Horizontal translations
        this->PeriodicIDs[1] = std::make_pair(-1, -1); // Vertical translations
    }

    virtual void setArea()
    {
        this->Area = boost::math::constants::pi<T>() * Radius * Radius;
    }
    virtual void setPerimeter()
    {
        this->Perimeter = static_cast<T>(2) * boost::math::constants::pi<T>() * Radius;
    }

    virtual Point<T> getPointInside()
    {
        //Random<T> Angle((T)0, (T)2 * boost::math::constants::pi<T>());
        //Random<T> Distance((T)0, (T)1);
        //T r = Radius * sqrt(Distance.getUniformRandom());
        //T a = Angle.getUniformRandom();
        //return Point<T>(r * cos(a), r * sin(a));
	return Point<T>(0, 0);
    }

    virtual bool isPointInside(const Point<T> &P)
    {
        return (P.x - Center.x) * (P.x - Center.x) + (P.y - Center.y) * (P.y - Center.y) < Radius * Radius;
    }

    // Currently edges are defined by single points, not curves so there is no implementation  for this
    // For complete generality, it is required to have an "identification" mechanism between arbitrary objects (cell edge is used for PBCs)
    virtual void setCellEdge()
    {
        this->CellEdge = std::make_pair(Point<T>(-1, -1), Point<T>(-1, -1));
    }

    virtual void setBoundingBox()
    {
        T val = Radius * tan(boost::math::constants::pi<T>() / (T)4);
        this->BoundingBox = std::make_pair(Point<T>(Center.x - val, Center.y - Radius), Point<T>(Center.x + val, Center.y + Radius));
    }

    virtual void setVertices()
    {
        this->Vertices = {Point<T>(NAN, NAN)}; 
    }

    virtual void Collide(Particle<T> &P, const int &OID)
    {

        Point<T> Normal(P.Pos.x - Center.x, P.Pos.y - Center.y);
        Normal = Normal / Norm(Normal); // Unit normal vector

        if (Dot<T>(P.Vel, Normal) < 0)
        {

            Point<T> V1 = P.Pos - Center;
            T V2 = Dot<T>(P.Vel, V1);
            T V3 = Dot<T>(V1, V1) - Radius * Radius;

            T Delta = V2 * V2 - V3;

            if (Delta > (T)0)
            {
                T TimeToCollision = -V2 - sqrt(Delta);
                P.UpdateState(TimeToCollision, OID, -1);
                //P.BirkhoffPos = Point<T>(getArcDistance(P.Pos, P.SegmentID), getIncidenceAngle(P.Vel, P.SegmentID));
            }
        }
    }

    virtual Point<T> getInwardNormal(const Particle<T>& P)
    {
        Point<T> Normal = P.Pos - Center;
        return Normal / Norm<T>(Normal);
    }

    virtual void Write(const std::string &Destination)
    {
    }

protected:
    T Radius;
    Point<T> Center;
};

#endif //CIRCLE_H_INCLUDED
