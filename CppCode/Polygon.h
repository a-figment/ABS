#ifndef POLYGON_H_INCLUDED
#define POLYGON_H_INCLUDED

#include <iomanip>
#include "Obstacle.h"
#include "LineSegment.h"
#include "Geometry.h"
#include "Particle.h"
#include "Random.h"

template <typename T>
class Polygon : public Obstacle<T>
{

public:
  Polygon() = default;
  Polygon(std::vector<T> xvertices, std::vector<T> yvertices, const std::vector<int> &pids, const bool &isInternal);
  //Polygon(std::vector<LineSegment<T>> segments, const std::vector<int> &pids);
  Polygon(std::vector<LineSegment<T>> segments, const std::vector<int> &pids);

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
  std::vector<Point<T>> getVertices() const { return xyVertices; }
  // std::pair<Point<T>, Point<T>> getBoundingBox() const { return BoundingBox; }

  virtual void Write(const std::string &Destination);
  bool isVertexCollision(const Point<T> &Pos);

private:
  int PIP(const std::vector<T> xVertices, const std::vector<T> yVertices, const Point<T> &P);
  // void setBoundingBox();
  void SegmentCollide(const LineSegment<T> &L, Particle<T> &P, T &ttc, const size_t &i, const int &OID);

  std::vector<int> PIDs;
  std::vector<T> xVertices;
  std::vector<T> yVertices;
  std::vector<Point<T>> xyVertices;
  std::vector<LineSegment<T>> Segments;
  // std::pair<Point<T>, Point<T>> BoundingBox;
  uint nVertices;
};

template <typename T>
Polygon<T>::Polygon(std::vector<T> xvertices, std::vector<T> yvertices, const std::vector<int> &pids, const bool &isinternal) : Obstacle<T>()
{
  assert(pids.size() % 2 == 0);
  assert(pids.size() != 0);
  assert(xvertices.size() == yvertices.size());
  assert(xvertices.size() > 1 || yvertices.size() > 1);

  xVertices = xvertices;
  yVertices = yvertices;
  PIDs = pids;

  this->isInternal = isinternal;
  nVertices = xVertices.size();
  this->nBoundaryComponents = nVertices;

  xyVertices.resize(nVertices);
  Segments.resize(nVertices);
  for (size_t i = 0; i < nVertices; ++i)
  {
    xyVertices[i] = Point<T>(xVertices[i], yVertices[i]);
    int j = (i + 1) % nVertices; // Assuming closed polygon
    Segments[i] = LineSegment<T>(Point<T>(xVertices[i], yVertices[i]), Point<T>(xVertices[j], yVertices[j]));
  }

  setPerimeter();
  setBoundingBox();
  setArea();
  setCumulativePerimeter();
  setPeriodicIDs(PIDs);
  setCellEdge();
  setVertices();
}

// @brief Here @Polygon is a collection of line segments, not all connected at their end points. 
// Assume only one side of each segment is accessible by particles.
template <typename T>
Polygon<T>::Polygon(std::vector<LineSegment<T>> segments, const std::vector<int> &pids) : Obstacle<T>()
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

  //xVertices.resize(nVertices);
  //yVertices.resize(nVertices);
  //xyVertices.resize(nVertices);

  for (size_t i = 0; i < Segments.size(); i++) {

    //xVertices[i]   = Segments[i].P.x;
    //xVertices[i+1] = Segments[i].Q.x;

    //yVertices[i]   = Segments[i].P.y;
    //yVertices[i+1] = Segments[i].Q.y;

    //xyVertices[i]   = Point<T>(Segments[i].P.x, Segments[i].P.y);
    //xyVertices[i+1] = Point<T>(Segments[i].Q.x, Segments[i].Q.y);

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
Point<T> Polygon<T>::getPointInside()
{
  Random<T> rx;
  Random<T> ry;

  Point<T> testPoint;
  bool isInside = false;
  while (!isInside)
  {
    testPoint = Point<T>(rx.getUniformRandom(this->BoundingBox.first.x, this->BoundingBox.second.x), 
			 ry.getUniformRandom(this->BoundingBox.first.y, this->BoundingBox.second.y));
    isInside = isPointInside(testPoint);
  }
  return testPoint;
}

// TODO: Implement general case of filling vector of pairs
template <typename T>
void Polygon<T>::setPeriodicIDs(const std::vector<int> &ids)
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
void Polygon<T>::setCellEdge()
{
  int lID = this->PeriodicIDs[0].second; // Left id
  int rID = this->PeriodicIDs[0].first;  // Right id

  int uID = this->PeriodicIDs[1].second; // Up id
  int dID = this->PeriodicIDs[1].first;  // Down id


  T xLeft  = Segments[lID].P.x;
  T xRight = Segments[rID].P.x;

  T xUp    = Segments[uID].P.x;
  T xDown  = Segments[dID].P.x;
  //std::cout << xLeft << " " << xRight << " " << xUp << " " << xDown << std::endl;

  //T xLeft  = xVertices[lID];
  //T xRight = xVertices[rID];

  //T xUp    = xVertices[uID];
  //T xDown  = xVertices[dID];

  this->CellEdge = std::make_pair(Point<T>(xLeft, xRight), Point<T>(xUp, xDown)); // Vertical translations not implemented
}

template <typename T>
void Polygon<T>::Collide(Particle<T> &P, const int &OID)
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

    const T VelocityOnNormal = Dot<T>(Q.Vel, L.InwardNormal);
    if (VelocityOnNormal > static_cast<T>(0)) // Particle<T> heading toward wall
      continue;

    const T v2 = Dot<T>(L.P - Q.Pos, L.InwardNormal);
    if (v2 > static_cast<T>(0)) // Particle<T> in front of wall
      continue;

    T TimeToCollision = v2 / VelocityOnNormal;
    Point<T> CollisionPoint(Q.Pos + TimeToCollision * Q.Vel);
    T DistanceToCollision = Distance2D<T>(CollisionPoint, L.MidPoint);
    if (DistanceToCollision > L.HalfLength)
      continue;

    if (TimeToCollision < ttc)
    { // Always true for first collision found

      Solution = P; // Re-assign - Solution may be updated more than once
      Solution.UpdateState(TimeToCollision, OID, i);
      Solution.BirkhoffPos = Point<T>(getArcDistance(P.Pos, P.SegmentID), getIncidenceAngle(P.Vel, P.SegmentID));
      ttc = TimeToCollision;
    }
  }
  P = Solution; // Guarantees no collisionFound 
}

template <typename T>
T Polygon<T>::getArcDistance(const Point<T> &P, const int &ID)
{
  return this->CumulativePerimeter[ID] + Distance2D<T>(xyVertices[ID], P);
}

template <typename T>
T Polygon<T>::getIncidenceAngle(const Point<T> &V, const int &ID)
{
  Point<T> Normal = Segments[ID].InwardNormal;
  Point<T> Velocity(-V.x, -V.y); // Incoming velocity
  return acos(Dot<T>(Normal, Velocity) / (Norm<T>(Normal) * Norm<T>(Velocity)));
}

template <typename T>
Point<T> Polygon<T>::getInwardNormal(const Particle<T>& P)
{
  return Segments[P.SegmentID].InwardNormal;
}

template <typename T>
void Polygon<T>::setBoundingBox()
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
void Polygon<T>::setPerimeter()
{
  this->Perimeter = static_cast<T>(0);
  for (size_t i = 0; i < nVertices; ++i)
    this->Perimeter += Distance2D<T>(xyVertices[i], xyVertices[(i + 1) % nVertices]);
}

template <typename T>
void Polygon<T>::setVertices() { this->Vertices = xyVertices; }

template <typename T>
void Polygon<T>::setArea()
{
  this->Area = static_cast<T>(0);
  for (size_t i = 0; i < nVertices; ++i)
    this->Area += ((xVertices[i] * yVertices[(i + 1) % nVertices]) - (xVertices[(i + 1) % nVertices] * yVertices[i]));
  this->Area *= (static_cast<T>(1) / static_cast<T>(2));
}

template <typename T>
void Polygon<T>::setCumulativePerimeter()
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
bool Polygon<T>::isPointInside(const Point<T> &P)
{
  return PIP(xVertices, yVertices, P) == 1;
}

template <typename T>
bool Polygon<T>::isVertexCollision(const Point<T> &Pos)
{
  for (const auto p : xyVertices)
  {
    if (eq_eps(p.x, Pos.x) && eq_eps(p.y, Pos.y))
      return true;
  }
  return false;
}

/*
 * @brief Point in polygon test
 * @return -1 Point OUTSIDE polygon
 *          0 Point ON polygon (edge/vertex)
 *          1 Point INSIDE polygon
 * @source https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html
 */
template <class T>
int Polygon<T>::PIP(const std::vector<T> xVertices, const std::vector<T> yVertices, const Point<T> &P)
{
  int i, j, c = 0;
  for (i = 0, j = static_cast<int>(xVertices.size()) - 1; i < static_cast<int>(xVertices.size()); j = i++)
    if (((yVertices[i] > P.y) != (yVertices[j] > P.y)) && (P.x < (xVertices[j] - xVertices[i]) * (P.y - yVertices[i]) / (yVertices[j] - yVertices[i]) + xVertices[i]))
      c = !c;
  return c;
}

// @brief Write polygon data (vertices) to file
    template <typename T>
void Polygon<T>::Write(const std::string &Destination)
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
std::ostream &operator<<(std::ostream &os, const Polygon<T> &P)
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
