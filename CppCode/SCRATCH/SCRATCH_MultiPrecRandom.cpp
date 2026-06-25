#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/random.hpp>

int main()
{
    using namespace boost::multiprecision;
    //    using namespace boost::random;
    //
    // Generate some distruted values:
    //
    boost::random::uniform_real_distribution<boost::multiprecision::cpp_dec_float_50> ur(-20, 20);
    boost::random::independent_bits_engine<boost::random::mt19937, std::numeric_limits<boost::multiprecision::cpp_dec_float_50>::digits, boost::multiprecision::cpp_int> gen;
    //
    // Generate some values:
    //
    std::cout << std::setprecision(50);
    for (unsigned i = 0; i < 20; ++i)
        std::cout << ur(gen) << std::endl;
    return 0;
}

// #include <boost/multiprecision/cpp_int.hpp>
// #include <boost/multiprecision/cpp_dec_float.hpp>
// #include <boost/random.hpp>

// typedef boost::multiprecision::cpp_dec_float<50> FLOAT;
// typedef boost::random::independent_bits_engine<boost::mt19937, std::numeric_limits<FLOAT>::digits, boost::multiprecision::cpp_int> gen;

// template <class T>
// inline T UniformRandom(gen &genny)
// {
//     // boost::random::random_device rd;
//     // boost::random::independent_bits_engine<boost::mt19937, std::numeric_limits<FLOAT>::digits, boost::multiprecision::cpp_int> gen;
//     // random_device rd;
//     // random::random_device rd;
//     // boost::random::mt19937 genny(gen());
//     // boost::random::random_number_generator<T> rng;
//     boost::random::uniform_real_distribution<T> dis;
//     return dis(genny);
// }

// int main()
// {
//     // fine
//     // double r1 = UniformRandom<double>(0.0L, 1.0L);
//     // not fine
//     // FLOAT r2;
//     // r2.
//     // FLOAT r2 = UniformRandom<FLOAT>(FLOAT(0.), FLOAT(1.));
//     std::cout << std::setprecision(50);
//     // gen genny(time(0));
//     gen genny;

//     for (unsigned i = 0; i < 3; ++i)
//     {
//         FLOAT p = UniformRandom<FLOAT>(genny);
//         // std::cout << UniformRandom<FLOAT>(genny) << std::endl;
//     }

//     return 0;
// }