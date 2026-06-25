#ifndef FLOATINGPOLYGON_H_INCLUDED
#define FLOATINGPOLYGON_H_INCLUDED

#include <iomanip>
#include "Obstacle.h"
#include "LineSegment.h"
#include "Geometry.h"
#include "Particle.h"
#include "Random.h"

template <typename T>
class FloatingPolygon : public Obstacle<T>
{

public:
  FloatingPolygon() = default;
  FloatingPolygon(std::vector<LineSegment<T>> segments, const std::vector<int> &pids);

  virtual void Collide(Particle<T> &P, const int &OID);
  virtual void setCumulativePerimeter();
  virtual void setPerimeter();
  virtual void setArea();
  virtual void setPeriodicIDs(const std::vector<int> &ids);
  virtual void setCellEdge();
  virtual void setBoundingBox();
  virtual void setVertices();

  virtual Point<T> getPointInside();
  virtual bool isPointInside(const Point<T> &P);

  virtual Point<T> getInwardNormal(const Particle<T>& P);
  virtual T getArcDistance(const Point<T> &P, const int &ID);
  virtual T getIncidenceAngle(const Point<T> &V, const int &ID);
  virtual void Write(const std::string &Destination);
  std::vector<Point<T>> getVertices() const { return xyVertices; }

private:
  std::vector<int> PIDs;
  std::vector<T> xVertices;
  std::vector<T> yVertices;
  std::vector<Point<T>> xyVertices;
  std::vector<LineSegment<T>> Segments;
  // std::pair<Point<T>, Point<T>> BoundingBox;
  uint nVertices;
};


// @brief Here @Polygon is a collection of line segments, not all connected at their end points. 
// Either side of the segment is accessible by particles
template <typename T>
FloatingPolygon<T>::FloatingPolygon(std::vector<LineSegment<T>> segments, const std::vector<int> &pids) : Obstacle<T>()
{
  assert(pids.size() % 2 == 0);
  assert(pids.size() != 0);
  //assert(xvertices.size() == yvertices.size());
  //assert(xvertices.size() > 1 || yvertices.size() > 1);
  assert(segments.size() >= 1);

  Segments = segments;
  PIDs = pids;
  this->isInternal = false;
  nVertices = 2*Segments.size(); // Two for each segment. Assuming each segment is not connected
			   
  //this->nBoundaryComponents = nVertices;
  this->nBoundaryComponents = segments.size();

  xVertices.reserve(nVertices);
  yVertices.reserve(nVertices);
  xyVertices.reserve(nVertices);

  for (size_t i = 0; i < Segments.size(); i++) {

    xVertices.push_back(Segments[i].P.x);
    xVertices.push_back(Segments[i].Q.x);

    yVertices.push_back(Segments[i].P.y);
    yVertices.push_back(Segments[i].Q.y);

    xyVertices.push_back(Point<T>(Segments[i].P.x, Segments[i].P.y));
    xyVertices.push_back(Point<T>(Segments[i].Q.x, Segments[i].Q.y));

  }

  setPerimeter();
  setBoundingBox();
  setArea();
  setCumulativePerimeter();
  setPeriodicIDs(PIDs);
  setCellEdge();
  setVertices();
}

// @brief Rejection samples the bounding box of the polygon with uniform x/y coordinates
template <typename T>
Point<T> FloatingPolygon<T>::getPointInside()
{
    return Point<T>(0,0);
}

// TODO: Implement general case of filling vector of pairs
template <typename T>
void FloatingPolygon<T>::setPeriodicIDs(const std::vector<int> &ids)
{
  // this->PeriodicIDs.resize(ids.size() / 2);
  this->PeriodicIDs.resize(4); // If only L/R pids are specified, fill the remainder by default
  if (ids.size() == 4)         // guaranteed by a precious assert
  {
    this->PeriodicIDs[0] = std::make_pair(ids[0], ids[1]); // Horizontal translations
    this->PeriodicIDs[1] = std::make_pair(ids[2], ids[3]); // Vertical translations
  }
  else
  {
    this->PeriodicIDs[0] = std::make_pair(ids[0], ids[1]);
    this->PeriodicIDs[1] = std::make_pair(-2, -2);
  }
}

template <typename T>
void FloatingPolygon<T>::setCellEdge()
{
    this->CellEdge = std::make_pair(Point<T>(0, 0), Point<T>(0, 0)); // Vertical translations not implemented
}

template <typename T>
void FloatingPolygon<T>::Collide(Particle<T> &P, const int &OID)
{
    const Particle<T> Q = P;  // Used as a temporary reference to the present state of P
    Particle<T> Solution = P; // Solution is updated from the state P for each time there is a valid event
    Solution.collisionFound = false; // Unsure if needed -- sanity.

    T ttc = INFINITY;
    for (size_t i = 0; i < Segments.size(); ++i)
    {
	// int j = (i + 1) % nVertices; // Assuming closed polygon
	// LineSegment<T> L(Point<T>(xVertices[i], yVertices[i]), Point<T>(xVertices[j], yVertices[j]));
	LineSegment<T> L = Segments[i];

	if (Q.SegmentID == static_cast<int>(i)) // Avoid glancing collisions
	    continue;

	Point<T> v0 = L.P - P.Pos; // Vector running from position to point P of L
	Point<T> v2 = L.Q - P.Pos; // Vector running from position to point Q of L
	bool passesThroughSegment = VectorBetweenVectors(v0, P.Vel, v2); 

	if (!passesThroughSegment)
	    continue;


	// Which normal vector faces toward the particle?
	LineSegment<T> L1(L.Q, L.P);

	T V = 0;
	T VelocityOnNormal = 0;
	if (Dot<T>(L.InwardNormal, Q.Vel) < 0) {
	    VelocityOnNormal = Dot<T>(Q.Vel, L.InwardNormal);
	    V = Dot<T>(L.P - Q.Pos, L.InwardNormal);

	    //std::cout << "Normal: " << L.InwardNormal << std::endl;
	} else {
	    VelocityOnNormal = Dot<T>(Q.Vel, L1.InwardNormal);
	    V = Dot<T>(L1.P - Q.Pos, L1.InwardNormal);
	    //std::cout << "Normal: " << L1.InwardNormal << std::endl;
	}


	T TimeToCollision = V / VelocityOnNormal;
	//std::cout << TimeToCollision << std::endl;
	Point<T> CollisionPoint(Q.Pos + TimeToCollision * Q.Vel);
	T DistanceToCollision = Distance2D<T>(CollisionPoint, L.MidPoint);
	if (DistanceToCollision > L.HalfLength)
	    continue;

	//std::cout << "TTC = " << TimeToCollision << std::endl;
	if (TimeToCollision < ttc)
	{ // Always true for first collision found

	    Solution = P; // Re-assign - Solution may be updated more than once
	    Solution.UpdateState(TimeToCollision, OID, i);
	    Solution.BirkhoffPos = Point<T>(-1, -1);
	    ttc = TimeToCollision;
	}
    }
    P = Solution; // Guarantees no collisionFound 
}

template <typename T>
T FloatingPolygon<T>::getArcDistance(const Point<T> &P, const int &ID)
{
    return this->CumulativePerimeter[ID] + Distance2D<T>(xyVertices[ID], P);
}

    template <typename T>
T FloatingPolygon<T>::getIncidenceAngle(const Point<T> &V, const int &ID)
{
    Point<T> Normal = Segments[ID].InwardNormal;
    Point<T> Velocity(-V.x, -V.y); // Incoming velocity
    return acos(Dot<T>(Normal, Velocity) / (Norm<T>(Normal) * Norm<T>(Velocity)));
}

    template <typename T>
Point<T> FloatingPolygon<T>::getInwardNormal(const Particle<T>& P)
{
    return Segments[P.SegmentID].InwardNormal;
}

    template <typename T>
void FloatingPolygon<T>::setBoundingBox()
{
    auto xMinMax = minmax_element(xVertices.begin(), xVertices.end());
    auto yMinMax = minmax_element(yVertices.begin(), yVertices.end());
    T xMax = *xMinMax.second;
    T yMax = *yMinMax.second;
    T xMin = *xMinMax.first;
    T yMin = *yMinMax.first;
    this->BoundingBox = std::make_pair(Point<T>(xMin, yMin), Point<T>(xMax, yMax));
}

    template <typename T>
void FloatingPolygon<T>::setPerimeter()
{
    this->Perimeter = static_cast<T>(0);
    for (size_t i = 0; i < nVertices; ++i)
	this->Perimeter += Distance2D<T>(xyVertices[i], xyVertices[(i + 1) % nVertices]);
}

template <typename T>
void FloatingPolygon<T>::setVertices() { this->Vertices = xyVertices; }

    template <typename T>
void FloatingPolygon<T>::setArea()
{
    this->Area = static_cast<T>(0);
    for (size_t i = 0; i < nVertices; ++i)
	this->Area += ((xVertices[i] * yVertices[(i + 1) % nVertices]) - (xVertices[(i + 1) % nVertices] * yVertices[i]));
    this->Area *= (static_cast<T>(1) / static_cast<T>(2));
}

    template <typename T>
void FloatingPolygon<T>::setCumulativePerimeter()
{
    this->CumulativePerimeter.reserve(nVertices);
    this->CumulativePerimeter[0] = static_cast<T>(0);
    for (size_t i = 0; i < nVertices - 1; ++i)
    {
	if (!((size_t)PIDs[0] == i || (size_t)PIDs[1] == i))
	    this->CumulativePerimeter.push_back(Distance2D<T>(xyVertices[i], xyVertices[i + 1])); // This gets the total perimeter
    }
    this->CumulativePerimeter = CumSum<T>(this->CumulativePerimeter);
}

    template <class T>
bool FloatingPolygon<T>::isPointInside(const Point<T> &P)
{
    //return PIP(xVertices, yVertices, P) == 1;
    return true;
}

//template <typename T>
//bool FloatingPolygon<T>::isVertexCollision(const Point<T> &Pos)
//{
//  for (const auto p : xyVertices)
//  {
//    if (eq_eps(p.x, Pos.x) && eq_eps(p.y, Pos.y))
//      return true;
//  }
//  return false;
//}

/*
 * @brief Point in polygon test
 * @return -1 Point OUTSIDE polygon
 *          0 Point ON polygon (edge/vertex)
 *          1 Point INSIDE polygon
 * @source https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html
 */
//template <class T>
//int FloatingPolygon<T>::PIP(const std::vector<T> xVertices, const std::vector<T> yVertices, const Point<T> &P)
//{
//  int i, j, c = 0;
//  for (i = 0, j = static_cast<int>(xVertices.size()) - 1; i < static_cast<int>(xVertices.size()); j = i++)
//    if (((yVertices[i] > P.y) != (yVertices[j] > P.y)) && (P.x < (xVertices[j] - xVertices[i]) * (P.y - yVertices[i]) / (yVertices[j] - yVertices[i]) + xVertices[i]))
//      c = !c;
//  return c;
//}
//
// @brief Write polygon data (vertices) to file
    template <typename T>
void FloatingPolygon<T>::Write(const std::string &Destination)
{
    std::ofstream Writer(Destination, std::ios::out);
    Writer << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;
    if (Writer.is_open())
	for (size_t i = 0; i < xVertices.size(); ++i)
	    Writer << xVertices[i] << " " << yVertices[i] << std::endl;
    else
    {
	std::cout << "Issue opening file: " << Destination << std::endl;
    }
    Writer.close();
}

    template <typename T>
std::ostream &operator<<(std::ostream &os, const FloatingPolygon<T> &P)
{
    os << "(    x,      y     ) = \n";
    for (const auto p : P.getVertices())
	os << p << ",\n";
    os << "Bounding box = {" << P.getBoundingBox().first << ", " << P.getBoundingBox().second << "}.";
    return os;
}

// DEBUG MULTIPRECISION
// std::cout << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;
// std::cout << "Velocity on normal    = " << VelocityOnNormal << std::endl;
// std::cout << "Velocity behind wall  = " << v2 << std::endl;
// std::cout << "Time to collision     = " << TimeToCollision << std::endl;
// std::cout << "Distance to collision = " << DistanceToCollision << std::endl;
// std::cout << "Midpoint              = " << L.MidPoint << std::endl;
// std::cout << "Half length           = " << L.Length / static_cast<T>(2) << std::endl
//           << std::endl;

#endif //POLYGON_H_INCLUDED
