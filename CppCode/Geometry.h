#ifndef GEOMETRY_H_INCLUDED
#define GEOMETRY_H_INCLUDED

#include <cmath>
#include <vector>
#include "Precision.h"
#include "Point.h"
#include "Particle.h"
#include "OnlineStatistics.h"

// 2D Distance
template <typename T>
inline T Distance2D(const Point<T> &P, const Point<T> &Q)
{
    T dx = Q.x - P.x;
    T dy = Q.y - P.y;
    return sqrt(dx * dx + dy * dy);
}

template <typename T>
inline T Dot(const Point<T> &P, const Point<T> &Q) { return P.x * Q.x + P.y * Q.y; }

// Length w.r.t origin
template <typename T>
inline T Norm(const Point<T> &P) { return sqrt(P.x * P.x + P.y * P.y); }

// Cumulative sum
// @source https://stackoverflow.com/questions/3288422/how-to-calculate-the-cumulative-sum-for-a-vector-of-doubles-in-c
template <typename T>
inline std::vector<T> CumSum(const std::vector<T> &V)
{
    std::vector<T> Result(V.size());
    partial_sum(V.begin(), V.end(), Result.begin());
    return Result;
}

template <typename T>
inline bool isPointOnLine(const Point<T> &P, const Point<T> &Q, const Point<T> &PointOnLine, const T &tol) // returns true if S is on line P, Q
{
    Point<T> X = P - Q;
    Point<T> Y = PointOnLine - Q;
    T Area = X.x * Y.y - X.y * Y.x;
    return fabs(Area) < tol;
}

/*
 * @brief Returns a (Point) vector which has been reflected about @Normal 
 * @detailed Vout = Vin - 2 <Vin, n> n, where <,> is the dot product
 * @Normal A vector (inward normal for billiards)
 * @Vin Incoming velocity vector
 */
template <typename T>
inline Point<T> ReflectSpecular(const Point<T> &Normal, const Point<T> &Vin)
{
    Point<T> Vout = Vin - (static_cast<T>(2) * Dot<T>(Vin, Normal) * Normal);
    //return Vout;
    return Vout / Norm(Vout);
}

template <typename T>
inline Point<T> getPointOnVerticalLine(const Point<T> &StartPt, const T &Slope, const T &Vertical)
{
    return Point<T>(Vertical, Slope * (Vertical - StartPt.x) + StartPt.y);
}

// @brief Computes the distance from Particle P1 (prior to escape) to L where it is assumed that
template <typename T>
inline T getDistanceToEscape(const Particle<T> &P1, const T &VerticalLine)
{
    Point<T> LiftedPos(P1.Pos.x + (T)P1.CellID.x, P1.Pos.y + (T)P1.CellID.y); // Prior to escaping

    T l = VerticalLine;
    if (P1.Vel.x < static_cast<T>(0))
        l = -VerticalLine;
    Point<T> EscapePoint = getPointOnVerticalLine<T>(LiftedPos, P1.Vel.y / P1.Vel.x, l);

    return Distance2D<T>(LiftedPos, EscapePoint);
}


// @A, @B [point] Start and end points of a line segment
// @x [vector] A series of coordinates {A.x, @x[0], @x[1], ..., @B.x}
// @y [vector] A series of coordinates {A.y, @y[0], @y[1], ..., @B.y}
// @return [vector] Perpendicular distances from the curve (@x, @y) to the segment AB
template <typename T>
std::vector<T> getPerpendicularDistances(const Point<T>& A, const Point<T>& B, const std::vector<T>& x, const std::vector<T>& y) {
    assert(x.size() == y.size());

    // Start and end points 
    T P1x = A.x;
    T P1y = A.y;
    T P2x = B.x;
    T P2y = B.y;

    // Denominator of perp distance eq
    T den = sqrt((P2x - P1x)*(P2x - P1x) + (P2y - P1y)*(P2y - P1y));
    
    std::vector<T> perp_dists(x.size(),0);
    for (size_t j = 0; j < x.size(); ++j) { 

        // Raw data points between A and B
        T x0 = x[j];
        T y0 = y[j];
    
        T num = std::abs((P2x - P1x) * (P1y - y0) - (P1x - x0)*(P2y - P1y));
        perp_dists[j] = num / den;
    }

    return perp_dists;
}


// @brief Gets the point on segment AB which has perp. line passing through C
// @return [point] 
template <typename T>
Point<T> getPointToLineNormal(const Point<T>& A, const Point<T>& B, const T& x, const T& y) {
    Point<T> C(x, y);
    T px = B.x - A.x;
    T py = B.y - A.y;
    T dAB = px*px + py*py;
    T u = ((C.x - A.x) * px + (C.y - A.y) * py) / dAB;
    return Point<T>(A.x + u*px, A.y + u*py);
}

//@brief 
template <typename T>
std::vector<Point<T>> getPointNormals(const Point<T>& A, const Point<T>& B, const std::vector<T>& x, const std::vector<T>& y) {
    assert(x.size() == y.size());
    
    std::vector<Point<T>> normal_comp(x.size(),Point<T>(0,0));
    for (size_t j = 0; j < x.size(); ++j) { 
	Point<T> P = getPointToLineNormal(A, B, x[j], y[j]);
	T x_comp = x[j] - P.x;
	T y_comp = y[j] - P.y;
	normal_comp[j] = Point<T>(x_comp, y_comp);

    }
    return normal_comp;
}

template <typename T>
T getMaxPerpendicularDistances(const Point<T>& A, const Point<T>& B, const std::vector<T>& x, const std::vector<T>& y) {
    assert(x.size() == y.size());

    // Start and end points 
    T P1x = A.x;
    T P1y = A.y;
    T P2x = B.x;
    T P2y = B.y;

    // Denominator of perp distance eq
    T den = sqrt((P2x - P1x)*(P2x - P1x) + (P2y - P1y)*(P2y - P1y));
    
    std::vector<T> perp_dists(x.size(),0);
    for (size_t j = 0; j < x.size(); ++j) { 

        // Raw data points between A and B
        T x0 = x[j];
        T y0 = y[j];
    
        T num = std::abs((P2x - P1x) * (P1y - y0) - (P1x - x0)*(P2y - P1y));
        perp_dists[j] = num / den;
    }

    T M = *std::max_element(perp_dists.begin(), perp_dists.end());
    return M;
}


template <typename T>
std::vector<T> StatSummary(const std::vector<T>& v) {

    OnlineStatistics<T> OS;
    OS.Clear();

    for (size_t i = 0; i < v.size(); ++i) {
	OS.Update(v[i]);
    }

    return { OS.getMean(), OS.getVariance() };

}

// Checks if a point is "left" or "right" of a line segment AB
// @return 0  point is on line
//	   +1 point is left of line
//	   -1 point is right of line
template <typename T>
int PointLine(const Point<T>& A, const Point<T>& B, const Point<T> C) {
    T val = ((B.x - A.x)*(C.y - A.y) - (B.y - A.y)*(C.x - A.x));
    return val / std::abs(val);
}

// Counts the number of points (x, y) left of the segment AB 
// @return P.x = num left points, P.y = num right points
template <typename T>
Point<int> numLeftRight(const Point<T>& A, const Point<T>& B, const std::vector<T>& x, const std::vector<T>& y) {

    int count_left = 0;
    int count_right = 0;
    for (size_t i = 0; i < x.size(); ++i) {
	int result = PointLine(A, B, Point<T>(x[i], y[i]));

	if (result > 0) {
	    count_left++;
	} else if (result < 0) {
	    count_right++;
	}
    }

    return Point<int>(count_left, count_right);
}

// @source https://stackoverflow.com/questions/243945/calculating-a-2d-vectors-cross-product
template <typename T>
T CrossProduct(const Point<T>& v1, const Point<T>& v2) 
{
    return (v1.x*v2.y - v1.y*v2.x);
}

// @brief Checks if v2 points 'between' v1 and v3
// @source https://stackoverflow.com/questions/13640931/how-to-determine-if-a-vector-is-between-two-other-vectors
template <typename T>
bool VectorBetweenVectors(const Point<T>& v1, const Point<T>& v2, Point<T>& v3) 
{
    T AxB = CrossProduct(v1, v2);
    T AxC = CrossProduct(v1, v3);
    T CxB = CrossProduct(v3, v2);
    T CxA = CrossProduct(v3, v1);
    return ((AxB * AxC >= (T)0) && (CxB * CxA >= (T)0));
}

#endif // GEOMETRY_H_INCLUDED
