#ifndef OBSTACLE_H_INCLUDED
#define OBSTACLE_H_INCLUDED

#include <iomanip>
#include <vector>
#include <cassert>
#include <fstream>
#include "Geometry.h"
#include "Point.h"
#include "Particle.h"

// @brief Any object which a particle collides with
template <typename T>
class Obstacle
{
public:
  Obstacle() = default;

  virtual void setArea() = 0;
  virtual void setPerimeter() = 0;
  virtual void setCumulativePerimeter() = 0;
  virtual void setPeriodicIDs(const std::vector<int> &ids) = 0;
  virtual void setCellEdge() = 0;
  virtual void setBoundingBox() = 0;
  virtual void setVertices() = 0;

  virtual void Collide(Particle<T> &P, const int &OID) = 0;

  virtual void Write(const std::string &Destination) = 0;

  virtual T getArcDistance(const Point<T> &P, const int &ID) = 0;
  virtual T getIncidenceAngle(const Point<T> &V, const int &ID) = 0;

  virtual bool isPointInside(const Point<T> &P) = 0;
  virtual Point<T> getPointInside() = 0;
  virtual Point<T> getInwardNormal(const Particle<T>& P) = 0;

  // GETTERS
  bool getLocation() const { return isInternal; }
  std::vector<T> getCumulativePerim() const { return CumulativePerimeter; }
  std::vector<Point<T>> getVertices() const { return Vertices; }
  T getArea() const { return Area; }
  T getPerimeter() const { return Perimeter; }
  uint getNumComponents() const { return nBoundaryComponents; }
  std::vector<std::pair<int, int>> getPeriodicIDs() const { return PeriodicIDs; }
  std::pair<Point<T>, Point<T>> getCellEdge() const { return CellEdge; }
  std::pair<Point<T>, Point<T>> getBoundingBox() const { return BoundingBox; }

protected:
  bool isInternal;                              // Obstacles are either 'inside' or 'outside'. Overlapping internal obstacles are not allowed, but overlapping internal/external are
  T Area, Perimeter;                            //
  uint nBoundaryComponents;                     // Polygons and piecewise objects
  std::vector<T> CumulativePerimeter;           //
  std::vector<Point<T>> Vertices;               //
  std::vector<std::pair<int, int>> PeriodicIDs; // Connects periodic boundaries
  std::pair<Point<T>, Point<T>> CellEdge;       // {(x1,x2), (y1,y2)} location of periodic boundary edge
  std::pair<Point<T>, Point<T>> BoundingBox;
};

template <typename T>
std::ostream &operator<<(std::ostream &os, const Obstacle<T> &O)
{
  os << std::fixed << std::setprecision(10)
     << "Area = " << O.getArea() << "\n"
     << "Perimiter = " << O.getPerimeter() << std::endl;
  os << "Cumul. Perimeter = ";
  for (const auto x : O.getCumulativePerim())
    os << x << " ";
  os << std::endl;

  return os;
}

#endif //OBSTACLE_H_INCLUDED
