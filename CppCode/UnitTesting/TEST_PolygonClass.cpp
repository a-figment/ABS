
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE PolygonTest
#include <boost/test/unit_test.hpp>
#include "../Polygon.h"





// @sources: (Tolerance) https://www.boost.org/doc/libs/1_60_0/libs/test/doc/html/boost_test/testing_tools/extended_comparison/floating_point.html
BOOST_AUTO_TEST_CASE(PolygonClass)
{
    // Area tests
    Polygon P({ZERO, ONE / 3, ONE / 3}, {ZERO, ZERO, ONE});

    std::cout << P;
    // FLOAT exp_area(1. / 2.);
    FLOAT exp_area = FLOAT((ONE / 2)) * FLOAT((ONE / 3));
    std::cout << std::setprecision(std::numeric_limits<FLOAT>::digits10) << std::showpoint;

    std::cout << exp_area << std::endl;
    std::cout << P.getArea() << std::endl;
    BOOST_CHECK_EQUAL(exp_area, P.getArea());

    // Perimeter tests
    // FLOAT exp_perim = FLOAT(2.) + sqrt(FLOAT(2.));
    FLOAT exp_perim = ONE + (ONE / 3) + sqrt(ONE + (ONE / 3) * (ONE / 3));
    BOOST_CHECK_EQUAL(exp_perim, P.getPerimeter());

    // Cumulative perimeter tests
    std::vector<FLOAT> exp_cumulperim = {0, 1, 2};
    BOOST_TEST(exp_cumulperim == P.getCumulativePerim());

    // Arc distance
    RNG rng;
    for (size_t i = 0; i < P.getVertices().size(); ++i)
    {
        FLOAT r = UniformRandom<FLOAT>(FLOAT(0), FLOAT(1), rng);
        Point<FLOAT> CP(P.getVertices()[i].x + r * (P.getVertices()[(i + 1) % P.getVertices().size()].x - P.getVertices()[i].x),
                        P.getVertices()[i].y + r * (P.getVertices()[(i + 1) % P.getVertices().size()].y - P.getVertices()[i].y));
        FLOAT exp_s = d(P.getVertices()[i], CP) + P.getCumulativePerim()[i];
        BOOST_CHECK_EQUAL(exp_s, P.ArcDistance(CP, i));
    }
}

BOOST_AUTO_TEST_CASE(RectangleClass)
{
    // Area tests
    Rectangle R0(ONE / 2, ONE / 3);
    FLOAT exp_area0 = FLOAT(2) * ((ONE / 2) + (ONE / 3));
    BOOST_CHECK_EQUAL(exp_area0, R0.getArea());

    Rectangle R1(ONE, ONE);
    FLOAT exp_area = ONE;
    BOOST_CHECK_EQUAL(exp_area, R1.getArea());

    // Perimeter tests
    Rectangle R2(FLOAT(2) / 10, FLOAT(5) / 10);
    FLOAT exp_perim = FLOAT(2) * (FLOAT(2) / 10) + FLOAT(2) * (FLOAT(5) / 10);
    BOOST_CHECK_EQUAL(exp_perim, R2.getPerimeter());

    // Cumulative perimeter tests
    Rectangle R3(Point<FLOAT>(1.0L, 1.0L), 1.0L, 1.0L);
    std::vector<FLOAT> exp_cumulperim = {0.0L, 1.0L, 2.0L, 3.0L};
    BOOST_TEST(exp_cumulperim == R3.getCumulativePerim());

    // Arc distance
    RNG rng;
    for (size_t i = 0; i < R3.getVertices().size(); ++i)
    {
        FLOAT r = UniformRandom<FLOAT>(0L, 1L, rng);
        Point<FLOAT> CP(R3.getVertices()[i].x + r * (R3.getVertices()[(i + 1) % R3.getVertices().size()].x - R3.getVertices()[i].x),
                        R3.getVertices()[i].y + r * (R3.getVertices()[(i + 1) % R3.getVertices().size()].y - R3.getVertices()[i].y));
        FLOAT exp_s = d(R3.getVertices()[i], CP) + R3.getCumulativePerim()[i];
        BOOST_CHECK_EQUAL(exp_s, R3.ArcDistance(CP, i));
    }
}
