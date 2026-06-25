
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE GeometryTest
#include <boost/test/unit_test.hpp>
// #include "../Geometry.h"
#include "../Random.h"
#include <math.h>

// @sources: (Tolerance) https://www.boost.org/doc/libs/1_60_0/libs/test/doc/html/boost_test/testing_tools/extended_comparison/floating_point.html
BOOST_AUTO_TEST_CASE(Dist2D)
{
    std::cout << std::setprecision(std::numeric_limits<FLOAT>::digits10) << std::showpoint;

    for (int i = 0; i < 2; ++i)
    {
        RNG genny;

        FLOAT A = URPointInInterval<FLOAT>(-FLOAT(1000000), FLOAT(1000000), genny);
        FLOAT B = URPointInInterval<FLOAT>(-FLOAT(1000000), FLOAT(1000000), genny);
        std::cout << A << " " << B << std::endl;

        // Point<FLOAT> P1 = URPointInRect(A, B, genny);
        // Point<FLOAT> P2 = URPointInRect(A, B, genny);

        // FLOAT exp_length = sqrt((P2 - P1).x * (P2 - P1).x + (P2 - P1).y * (P2 - P1).y);

        // std::cout << P1 << P2 << std::endl;

        // BOOST_TEST(Distance2D(P1, P2), exp_length);
    }
}
