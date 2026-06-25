#ifndef PRECISION_INCLUDED
#define PRECISION_INCLUDED

// #include <mpfr.h>
// #include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/math/constants/constants.hpp>

namespace mp = boost::multiprecision;
typedef double float_dub;
typedef long double float_ldub;
typedef long long int int_ll;
typedef long int cid_int;

// typedef mp::cpp_dec_float<300> mp_backend;
// typedef mp::number<mp_backend, mp::et_off> FLOAT; // Expression templates off

typedef mp::number<mp::cpp_dec_float<50>, mp::et_off> float_50; 
//typedef mp::number<mp::cpp_int::backend_type, mp::et_off> int_ll;

typedef mp::number<mp::cpp_dec_float<100>, mp::et_off> float_100;
typedef mp::number<mp::cpp_dec_float<200>, mp::et_off> float_200;
typedef mp::number<mp::cpp_dec_float<300>, mp::et_off> float_300;
typedef mp::number<mp::cpp_dec_float<400>, mp::et_off> float_400;
typedef mp::number<mp::cpp_dec_float<500>, mp::et_off> float_500;
typedef mp::number<mp::cpp_dec_float<600>, mp::et_off> float_600;
typedef mp::number<mp::cpp_dec_float<700>, mp::et_off> float_700;
typedef mp::number<mp::cpp_dec_float<800>, mp::et_off> float_800;
typedef mp::number<mp::cpp_dec_float<900>, mp::et_off> float_900;
typedef mp::number<mp::cpp_dec_float<1000>, mp::et_off> float_1000;

// const FLOAT ONE = static_cast<FLOAT>(1);
// const FLOAT ZERO = static_cast<FLOAT>(0);

template <typename T>
bool eq_eps(const T &a, const T &b)
{
        return (fabs(a - b) < std::numeric_limits<T>::epsilon());
}

template <typename T>
bool gt_eps(const T &a, const T &b)
{
        return (a - b) > ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * std::numeric_limits<T>::epsilon());
}

template <typename T>
bool lt_eps(const T &a, const T &b)
{
        return (b - a) > ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * std::numeric_limits<T>::epsilon());
}

// bool definitelyGreaterThan(float a, float b, float epsilon)
// {
//     return (a - b) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
// }

// bool definitelyLessThan(float a, float b, float epsilon)
// {
//     return (b - a) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
// }

#endif // PRECISION_INCLUDED
