#ifndef SCRATCH_UNITTESTBOOST_INCLUDED
#define SCRATCH_UNITTESTBOOST_INCLUDED

#include <boost/multiprecision/cpp_dec_float.hpp>
// typedef double FLOAT;
typedef boost::multiprecision::cpp_dec_float_50 FLOAT;
const FLOAT ONE = FLOAT(1);

struct Rect
{
    Rect(const FLOAT &width, const FLOAT &height) : Width(width), Height(height){};
    FLOAT getArea() const
    {
        return FLOAT(Width * Height);
    }
    FLOAT Width, Height;
};

#endif
