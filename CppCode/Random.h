#ifndef RANDOM_H_INCLUDED
#define RANDOM_H_INCLUDED

#include "Geometry.h"
#include <random>
#include <boost/random.hpp>

// typedef boost::random::independent_bits_engine<boost::mt19937, std::numeric_limits<FLOAT>::digits, boost::multiprecision::number<boost::multiprecision::cpp_int::backend_type, boost::multiprecision::et_off>> generator;
// independent_bits_engine<boost::mt19937, std::numeric_limits<cpp_dec_float_50>::digits, cpp_int> gen(time(NULL));

// typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<50>, boost::multiprecision::et_off> float_50;   // Expression templates off
// typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<100>, boost::multiprecision::et_off> float_100;
// typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<200>, boost::multiprecision::et_off> float_200;
// typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<400>, boost::multiprecision::et_off> float_400;
// typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<800>, boost::multiprecision::et_off> float_800;

template <typename T>
class Random
{
public:
    using genny = boost::random::independent_bits_engine<boost::mt19937, std::numeric_limits<T>::digits, boost::multiprecision::cpp_int>;



    Random()
    {
        std::random_device rd;
        gen = genny(rd());
	//std::cout << "Seed fixed!" << std::endl;
        //gen = genny(1);
    }

    //Random(const T &a, const T &b)
    //{
    //    ur = boost::random::uniform_real_distribution<T>(a, b);
    //    std::random_device rd;
    //    gen = genny(rd());
    //    // gen = genny(rand());
    //    ur(gen);
    //}

    T getUniformRandom(const T &a, const T &b)
    {
        ur = boost::random::uniform_real_distribution<T>(a, b);
        return ur(gen);
    }

    T getNormalRandom(const T &mean, const T &var)
    {
        norm = boost::random::normal_distribution<T>(mean, var);
        return norm(gen);
    }

private:
    boost::random::uniform_real_distribution<T> ur;
    boost::random::normal_distribution<T> norm;
    genny gen;
};

// generator &gen()
// {
//     thread_local static generator genny(time(NULL));
//     return genny;
// }

// // @source https://stackoverflow.com/questions/686353/random-float-number-generation
// template <class T>
// inline T SimpleRandom(const T &lb, const T &ub)
// {
//     return lb + (rand()) / ((RAND_MAX / (ub - lb)));
// }

// // @brief Draws a number from [@lb, @ub] with uniform probability
// template <class T>
// inline T URPointInInterval(const T &lb, const T &ub, generator &Genny)
// {
//     boost::random::uniform_real_distribution<T> UNIF(lb, ub);
//     return UNIF(Genny);
// }

// // @brief Draws a point from a rectangle
// // @lc @uc lower and upper corners
// template <typename T> ////, class generator
// Point<T> URPointInRect(const Point<T> &lc, const Point<T> &uc)
// {
//     // generator Genny(time(NULL));
//     boost::random::uniform_real_distribution<T> xUNIF(lc.x, uc.x);
//     boost::random::uniform_real_distribution<T> yUNIF(lc.y, uc.y);
//     return Point<T>(xUNIF(gen()), yUNIF(gen()));
// }

// template <typename T> ////, class generator
// Point<T> URPointOnCircle()
// {
//     boost::random::uniform_real_distribution<T> UR(FLOAT(0), boost::math::constants::pi<T>());
//     T Theta = UR(gen());
//     return Point<T>(cos(Theta), sin(Theta));
// }

// @brief
// template <typename T>
// inline Point<T> URPointInRect(const T &a, const T &b, generator &Genny)
// {
//     boost::random::uniform_real_distribution<T> xUNIF(-a, a);
//     boost::random::uniform_real_distribution<T> yUNIF(-b, b);
//     return Point<T>(xUNIF(Genny), yUNIF(Genny));
// }

#endif // RANDOM_H_INCLUDED
