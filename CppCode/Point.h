#ifndef POINT_H_INCLUDED
#define POINT_H_INCLUDED
#include <iostream>
#include <iomanip>
#include <limits>
#include <vector>

template <typename T>
struct Point
{
    T x, y;
    Point<T>() = default;
    Point<T>(const T &X, const T &Y) : x(X), y(Y){};
};

template <typename T>
inline Point<T> operator+(const Point<T> &P, const Point<T> &Q) { return {P.x + Q.x, P.y + Q.y}; }

template <typename T>
inline Point<T> operator-(const Point<T> &P, const Point<T> &Q) { return {P.x - Q.x, P.y - Q.y}; }

template <typename T>
inline Point<T> operator/(const Point<T> &P, const T &a) { return {P.x / a, P.y / a}; }

template <typename T>
inline Point<T> operator*(const T &a, const Point<T> &P) { return {a * P.x, a * P.y}; }

template <typename T>
inline Point<T> operator*(const Point<T> &P, const Point<T> &Q) { return {P.x * Q.x, P.y * Q.y}; }

template <typename T>
std::ostream &operator<<(std::ostream &os, const Point<T> &P)
{
    return os << "\033[1;31m(" << P.x << ", " << P.y << ")\033[0m";
}
#endif // POINT_H_INCLUDED
