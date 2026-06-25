#ifndef ARCSEGMENT_H_INCLUDED
#define ARCSEGMENT_H_INCLUDED

#include <iomanip>
#include <vector>
#include "Geometry.h"
#include "Particle.h"
#include "Obstacle.h"

// @brief This class defines and 'arc segment' obstacle which is defined by a centre point
//	  (the centre of a full circle), radius, and two points defined along the perimeter of
//	  a circle. A collision event occurs when the velocity lies between the two vectors 
//	  drawn from the two arc segment points.
//	  We treat circular arc segments as a semi-circle, but with different conditions for when
//	  a collision occurs. 
//	  The segments are assumed to be symmatry about chord mid point (arc length equal either side)
// @WARNING. It is assumed that this obstacled is impacted from the 'inside'. For external scattering see Circle.h
template <typename T>
class ArcSegment : public Obstacle<T>
{

public:
    // ArcSegment(const T& radius, const Point<T>& center) : Radius(radius), Center(center) {}
    ArcSegment(const Point<T> &center, const Point<T>& arcpoint1, const Point<T>& arcpoint2) : Obstacle<T>()
    {
        Center = center; // CIRCLE CENTRE! Not the chord midpoint

	// AP1 and AP2 must be ordered appropriately! e.g. counter clockwise
	AP1 = arcpoint1;
	AP2 = arcpoint2;

	Radius = Distance2D<T>(AP1, Center);

	// Obstacle always hit from the inside. Direction is given by the normal 
	Midpoint = (AP1 + AP2) / (T)2;

	// Faces away from the arc (since we only collide from the inside)
        Direction = Point<T>(-(AP2.y - AP1.y), (AP2.x - AP1.x));
	Direction = Direction / Norm<T>(Direction);

	//std::cout << "Radius = " << Radius << std::endl;
	//std::cout << "Direction = " << Direction << std::endl;
	//std::cout << "Chord midpoint = " << Midpoint << std::endl << std::endl;

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

    // Not used (meaningless funcion)
    virtual void setArea()
    {
        this->Area = boost::math::constants::pi<T>() * Radius * Radius;
    }

    // Not implemented 
    virtual void setPerimeter()
    {
        this->Perimeter = static_cast<T>(2) * boost::math::constants::pi<T>() * Radius;
    }

    // Not used (meaningless funcion), not implemented
    virtual Point<T> getPointInside()
    {
        return Point<T>(0, 0);
    }

    // Used when generating initial conditions
    virtual bool isPointInside(const Point<T> &P)
    {
	bool awayFromDirection = (Dot<T>(Direction, P - Center) < 0);
	bool outsideSector = Distance2D<T>(P, Center) < Radius;
        return (awayFromDirection && outsideSector) || (Dot<T>(P - Center, Direction) > 0);
    }

    // Currently edges are defined by single points, not curves so there is no implementation  for this
    // For complete generality, it is required to have an "identification" mechanism between arbitrary objects (cell edge is used for PBCs)
    virtual void setCellEdge()
    {
        this->CellEdge = std::make_pair(Point<T>(-1, -1), Point<T>(-1, -1));
    }

    virtual void setBoundingBox()
    {
	//auto xMinMax = minmax_element();
	//auto yMinMax = minmax_element(yVertices.begin(), yVertices.end());
	//T xMax = *xMinMax.second;
	//T yMax = *yMinMax.second;
	//T xMin = *xMinMax.first;
	//T yMin = *yMinMax.first;
	//T minx = 
	//    if (Direction.y > 0) { // ASSUMING THERE ARE ONLY 'CUPS (U)' AND 'CAPS' AS ARCS

	//	T yMin = Center.y - Radis;
	//	T yMax = yMin; 

	//	T xMin = std::min(AP1.x, AP2.x);
	//	T xMax = std::max(AP1.x, AP2.x);

	//    }

	//}
        this->BoundingBox = std::make_pair(Point<T>(AP2.x, Midpoint.y), Point<T>(AP1.x, Midpoint.y + Radius));
    }

    virtual void setVertices()
    {
        this->Vertices = {Point<T>(NAN, NAN)}; 
    }

    virtual void Collide(Particle<T> &P, const int &OID)
    {

	// Particle currently outside the circle sector
	Point<T> PosToMidpoint = P.Pos - Midpoint;
	bool outsideArc = (Dot<T>(PosToMidpoint,  Direction) > 0); 

	// The following if-clauses are to check whether the particle will collide with the arc segment
	// If it does, we pretend that 
	if (outsideArc) { // Entering vicinity of arc from outside chord

	    // Particle velocity is between the vectors coming from the arc points
	    Point<T> v0 = AP1 - P.Pos; // Vector running from position to ap1
	    Point<T> v2 = AP2 - P.Pos; // Vector running from position to ap2
	    bool willCollideFromOutside = VectorBetweenVectors(v0, P.Vel, v2); 
	    if (willCollideFromOutside) {
		//std::cout << "outside chord" << std::endl;

		// Pretend the collision is with a semi-circle
		Point<T> PosToCenter = P.Pos - Center;
		T B = Dot<T>(P.Vel, PosToCenter);
		T C = Dot<T>(PosToCenter, PosToCenter) - Radius * Radius;
                T TimeToCollision = -B + sqrt(B*B - C);
                P.UpdateState(TimeToCollision, OID, -1);
		//std::cout << P << std::endl;
	    }

	} else { // Entering from within the sector: P.Pos.y - C.y > 0 for upper segments
	    
	    //Point<T> v0 = P.Pos - AP1; // Vector running from ap1 to position
	    //Point<T> v2 = P.Pos - AP2; // Vector running from ap2 to position
	    Point<T> v0 = AP1 - P.Pos; // Vector running from position to ap1
	    Point<T> v2 = AP2 - P.Pos; // Vector running from position to ap2
	    bool willCollideFromOutside = VectorBetweenVectors(v0, P.Vel, v2); 
	    if (!willCollideFromOutside) { // Moving from point on arc to another point on arc
		// Pretend the collision is with a semi-circle
		Point<T> PosToCenter = P.Pos - Center;
		T B = Dot<T>(P.Vel, PosToCenter);
		T C = Dot<T>(PosToCenter, PosToCenter) - Radius * Radius;
                T TimeToCollision = -B + sqrt(B*B - C);
                P.UpdateState(TimeToCollision, OID, -1);
	    }
	}

    }

    virtual Point<T> getInwardNormal(const Particle<T>& P)
    {
        //Point<T> Normal = P.Pos - Center; // Collision happens outside the circle
        Point<T> Normal = Center - P.Pos;  // Collision from insice
        return Normal / Norm<T>(Normal);
    }

    virtual void Write(const std::string &Destination)
    {
    }

protected:
    T Radius;
    Point<T> Center;
    Point<T> Direction; 
    Point<T> Midpoint; 
    Point<T> AP1;
    Point<T> AP2;
};

#endif //CIRCLE_H_INCLUDED
