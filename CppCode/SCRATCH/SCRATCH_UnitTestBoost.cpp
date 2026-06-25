
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE RectTest
#include <boost/test/unit_test.hpp>
#include "SCRATCH_UnitTestBoost.h"
namespace utf = boost::unit_test;

// Failing
BOOST_AUTO_TEST_CASE(AreaTest1)
{
    // Area tests
    Rect R(ONE / 2, ONE / 3);
    FLOAT expected_area = (ONE / 2) * (ONE / 3);

    // std::cout << std::setprecision(std::numeric_limits<FLOAT>::digits10); // << std::showpoint;
    // std::cout << "Expected: " << expected_area << std::endl;
    // std::cout << "Actual  : " << R.getArea() << std::endl;

    BOOST_CHECK_EQUAL(expected_area, R.getArea());
    // BOOST_TEST(expected_area == R.getArea());
}

// Tolerance has no effect?
BOOST_AUTO_TEST_CASE(AreaTestTol, *utf::tolerance(1e-40))
{
    // Area tests
    Rect R(ONE / 2, ONE / 3);
    FLOAT expected_area = (ONE / 2) * (ONE / 3);

    std::cout << std::setprecision(std::numeric_limits<FLOAT>::digits10) << std::showpoint;
    std::cout << "Expected: " << expected_area << std::endl;
    std::cout << "Actual  : " << R.getArea() << std::endl;

    // BOOST_CHECK_EQUAL(expected_area, R.getArea());
    BOOST_TEST(expected_area == R.getArea());
}

// Passing
BOOST_AUTO_TEST_CASE(AreaTest2)
{
    Rect R(ONE / 7, ONE);
    FLOAT expected_area = (ONE / 7) * (ONE);
    BOOST_CHECK_EQUAL(expected_area, R.getArea());
}

// Passing
BOOST_AUTO_TEST_CASE(AreaTest3)
{
    Rect R(ONE, ONE);
    FLOAT expected_area = ONE;

    std::cout << std::setprecision(std::numeric_limits<FLOAT>::digits10) << std::showpoint;
    std::cout << "Expected: " << expected_area << std::endl;
    std::cout << "Actual  : " << R.getArea() << std::endl;

    BOOST_CHECK_EQUAL(expected_area, R.getArea());
}
