
#define BOOST_TEST_DYN_LINK        // this is optional
#define BOOST_TEST_MODULE RectangleTest
#include <boost/test/unit_test.hpp>
#include "../Polygon.h" // project being tested

BOOST_AUTO_TEST_CASE(my_boost_test)
{

    Rectangle R(1,1);

    FLOAT exp_area = 1.;

    BOOST_CHECK(exp_area == R.getArea());
}

