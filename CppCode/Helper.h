#ifndef HELPER_INCLUDED
#define HELPER_INCLUDED

#include "Ensemble.h"
#include <algorithm>

// ASSORTED HELPER FUNCTIONS WITH NO HOME
//
//
//
//

//#include <random>
//#include <boost/random.hpp>
//
//template <typename T> 
//inline std::vector<T> getNormalDistribution(const size_t& n, const T& mean, const T& var) {
//
//    boost::random::normal_distribution<T> norm = boost::random::normal_distribution<T>(mean, var);
//    std::random_device rd;
//    boost::random::independent_bits_engine<boost::mt19937, std::numeric_limits<T>::digits, boost::multiprecision::cpp_int> gen = boost::random::independent_bits_engine<boost::mt19937, std::numeric_limits<T>::digits, boost::multiprecision::cpp_int>(rd());
//
//    std::vector<T> result(n,0);
//    for (size_t i = 0; i < n; ++i) {
//	result[i] = norm(gen);
//	std::cout << result[i] << std::endl;
//    }
//    return norm(gen);
//
//}



// @brief Is value contained in @vect. Only for integer types!
template <typename T> 
inline bool In(const T& value, const std::vector<T>& vect) {
    return (std::find(vect.begin(), vect.end(), value) != vect.end());
}

// @brief Sequence of @length numbers on [@a, @b]. 
template <typename T>
inline std::vector<T> seq(const T &a, const T &b, const int &length)
{
    std::vector<T> result(length);
    for (int i = 0; i <= length - 1; ++i)
        result[i] = a + (b - a) / (length - static_cast<T>(1)) * static_cast<T>(i);
    return result;
}

// @brief Same as seq() above but with print option
template <typename T>
inline std::vector<T> pseq(const T &a, const T &b, const int &length)
{
    std::vector<T> result = seq<T>(a, b, length);
    std::cout << "seq(" << a << ", " << b << ", " << length << ") = ";
    for (int i = 0; i <= length - 1; ++i)
        std::cout << result[i] << " ";
    std::cout << std::endl;
    return result;
}

// =========================== OVERLOADS START ===========================
// @brief WRITE vector to a file
template <typename T>
std::ofstream &operator<<(std::ofstream &ofs, const std::vector<T> &v)
{
    ofs << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;
    for (size_t i = 0; i < v.size(); ++i)
        ofs << v[i] << " ";
    return ofs; 
}

// @brief PRINT vector
template <typename T>
std::ostream &operator<<(std::ostream &ofs, const std::vector<T> &v)
{
    for (size_t i = 0; i < v.size()-1; ++i)
        ofs << v[i] << ", ";
    ofs << v.back() << ".\n";
    return ofs; 
}


// @brief PRINT a vector of points. See @Point<> for << overload.
template <typename T>
std::ostream &operator<<(std::ostream &ofs, const std::vector<Point<T>> &v)
{
    for (size_t i = 0; i < v.size(); ++i)
        ofs << v[i] << " ";
    return ofs; 
}

// @brief WRITE vector of points (@x, @y) as a column
template <typename T>
std::ofstream &operator<<(std::ofstream &ofs, const std::vector<Point<T>> &v)
{
    for (size_t i = 0; i < v.size(); ++i)
        ofs << v[i].x << " " << v[i].y << "\n";
    return ofs; 
}

// =========================== OVERLOADS END ===========================




// @brief Gets @n rational multiples of pi either side of pi/2 not including pi/2
// e.g. n = 3 produces pi/5, pi/4, pi/3,< pi/2, > 2pi/3, 3pi/4, 4pi/5
template <typename T>
inline std::vector<T> getRationalMultipleOfPi(const uint& n) { 

    T PI = boost::math::constants::pi<T>();

    std::vector<uint> l_denom = seq<uint>(2, 2 + n, n + 1); // includes pi/2
    std::reverse(l_denom.begin(), l_denom.end());
    std::vector<uint> r_numer = seq<uint>(2, n+1, n);
    std::vector<uint> r_denom = seq<uint>(3, 2 + n, n);

    std::vector<T> result(2*n + 1, 0);
    for (size_t i = 0; i < result.size(); ++i) {
	if (i <= n) 
	    result[i] = (1/(T)l_denom[i]) * PI;
	else
	    result[i] = ((T)r_numer[i-n-1] / (T)r_denom[i-n-1]) * PI;
    }

    std::cout << "Angles in study: " << std::endl;
    for (size_t i = 0; i < l_denom.size(); ++i)
	std::cout << i << ": PI/" << l_denom[i] << ", ";
    std::cout << std::endl;
    for (size_t i = 0; i < r_numer.size(); ++i)
	std::cout << i + l_denom.size() << ": " << r_numer[i] << "PI/" << r_denom[i] << ", ";
    std::cout << std::endl;

    return result;
}

template <typename T>
inline T cot(const T& a) { return (cos(a) / sin(a)); }

template <typename T>
inline T IntegerPart(const T& v) { 
    if (v < 0)
	return std::ceil(v);
    else
	return std::floor(v);
}

// @brief All elements @v are greater than @value
template <typename T>
inline bool all_gt(const std::vector<T>& v, const T& value) { 
    for (size_t i = 0; i < v.size(); ++i) 
	if (v[i] < value) // One of the elements is less than (and not greater than) @value
	    return false;
    return true;
}

// @brief All elements @v are less than @value
template <typename T>
inline bool all_lt(const std::vector<T>& v, const T& value) { 
    for (size_t i = 0; i < v.size(); ++i) 
	if (v[i] > value) // One of the elements is greater than @value
	    return false;
    return true;
}

template <typename T> 
std::vector<T> Concat(const std::vector<T>& A, const std::vector<T>& B) {
    std::vector<T> AB = A;
    AB.insert(AB.end(), B.begin(), B.end());
    return AB;
}


// Is value in epsilon n.h about mid
template <typename T>
inline bool in_nh(const T& v, const T& mid, const T& eps) {
    return (v < mid + eps && v > mid - eps);
}

template <typename T>
inline bool in_nh_ublb(const T& v, const T& lb, const T& ub) {
    return (v < ub && v > lb);
}


// @brief Indexes @vect at @ind returning @vect[@ind]
template <typename T>
inline std::vector<T> index(const std::vector<T>& vect, const std::vector<size_t>& ind) {
    std::vector<T> result(ind.size());
    for (size_t i = 0; i < ind.size(); ++i)
	result[i] = vect[ind[i]];
    return result;

}

template <typename T>
inline std::vector<T> index(const std::vector<T>& vect, const size_t& from, const size_t& to) {
    //assert(to > from);
    std::vector<T> result(to - from);
    for (size_t i = from; i < to; ++i)
	result[i] = vect[i];
    return result;

}


// @brief Given a @vect of ascending values, find the index corresponding to @valuesToFind
template <typename T>
inline std::vector<size_t> getIndicesOfVectorValue(const std::vector<T> valuesToFind, const std::vector<T>& vect) {

    std::vector<size_t> Indices(valuesToFind.size(), 0);
    for (size_t i = 0; i < valuesToFind.size(); ++i) {
        for (size_t j = 0; j < vect.size(); ++j) {
            if ((vect[j]) >= valuesToFind[i]) {
                Indices[i] = j;
                break;
            }
        }
    }

   std::cout << "Indices found: ";
   for (size_t i = 0; i < Indices.size(); ++i)
       std::cout << Indices[i] << " ";
   std::cout << std::endl;
   return Indices; 

}

// @brief Degrees to radirm
template <typename T>
inline T toRad(const T &ang) { return ang * (boost::math::constants::pi<T>() / ((T)180)); }

// @brief Radirm to degrees
template <typename T>
inline T toDeg(const T &ang) { return ang * (((T)180) / boost::math::constants::pi<T>()); }

// @brief Vector of @nSteps running from {1,...,@totalTime} in exponential form
//	  Equivalent to: exp(seq(0, log(@totalTime), length=@nSteps))
template <typename T>
inline std::vector<T> getExponentialTimes(const T &totalTime, const size_t &nSteps)
{
    std::cout << "Using log time with " << nSteps << " steps." << std::endl;
    T step = log(totalTime);
    std::vector<T> Time(nSteps, 0);
    for (size_t i = 0; i < nSteps; ++i)
    {
        T val = (step / (nSteps - static_cast<T>(1))) * static_cast<T>(i);
        Time[i] = exp(val) - 1;
    }
    return Time;
}

// @brief Vector of @nSteps running from {1,...,@totalTime} in exponential form
//	  Equivalent to: exp(seq(0, log(@totalTime), length=@nSteps))
template <typename T>
inline std::vector<T> getLinearTimes(const T &totalTime, const T &dt)
{
    std::cout << "Using linear times." << std::endl;
    std::vector<T> Time(static_cast<size_t>(totalTime)); // Continuous times to observe particle
    for (size_t i = 1; i <= Time.size(); ++i)
        Time[i - 1] = dt * static_cast<T>(i);
    return Time;
}


template <typename T>
inline std::vector<T> multiply(const std::vector<T> &v, const T &a)
{
    std::vector<T> result(v.size());
    for (size_t i = 0; i < v.size(); ++i)
        result[i] = v[i] * a;
    return result;
}

// @brief Updates @Ens with the displacements of a single @Traj
// @MomentStep Moments of @Ens computed as @MomentStep, @MomentStep^2,... @MomentStep^@nMoments
template <typename T>
inline void FillMoments(std::vector<Ensemble<T>> &Ens, const Trajectory<T> &Traj, const T &MomentStep, const int &nMoments)
{
    std::vector<std::vector<T>> DXp = Traj.getDisplacements_x(MomentStep, nMoments);
    for (int j = 0; j < nMoments; ++j)
        Ens[j].Push(DXp[j]);
}

//template <typename T>
//void WriteVector(const std::vector<T> &v, std::ofstream& Writer)
//{
//    Writer << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;
//
//    if (Writer.is_open()) {
//        for (size_t i = 0; i < v.size(); ++i)
//            Writer << v[i] << " ";
//        Writer << std::endl;
//    } else {
//        std::cout << "Issue opening file: " << std::endl;
//    }
//
//}


template <typename T>
void WritePositionVector(const std::string &Destination, const std::vector<Point<T>> &pv)
{
    std::ofstream Writer(Destination, std::ios::out | std::ios::app);
    Writer << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;
    if (Writer.is_open())
    {
        for (size_t i = 0; i < pv.size(); ++i)
            Writer << pv[i].x << " " << pv[i].y << '\n';
    }
    else
    {
        std::cout << "Issue opening file: " << Destination << std::endl;
    }
    Writer.close();
}

// @brief Takes "x" or "y" (@coord) component of a collection of @Points
template <typename T>
inline std::vector<T> getCoordFromVectorPoint(const std::vector<Point<T>> Points, const std::string &coord)
{
    std::vector<T> coordValue(Points.size());
    if (coord == "x")
        for (size_t i = 0; i < Points.size(); ++i)
            coordValue[i] = Points[i].x;
    else if (coord == "y")
        for (size_t i = 0; i < Points.size(); ++i)
            coordValue[i] = Points[i].y;
    else
        std::cout << "ERR: Misc::getCoordFromVectorPoint - Invalid 'coord' string." << std::endl;
    return coordValue;
}


// ========================================= CHANNEL CREATORS START =========================================
// Various functions which create a @Polygon object intended for representing various polygonal 
// billiards. Definitions are self-explanatory based on parameters accepted by each function.
//

template <typename T>
inline Polygon<T> getNPChannel_ByAngle(const T &angle_top, const T &angle_bot, const T &dx, const T &d)
{
    std::cout << "Creating non-parallel channel with bottom angle = " << toDeg(angle_bot) << " and top angle = " << toDeg(angle_top) << std::endl;

    T dy_b = dx * (cos(angle_bot / static_cast<T>(2)) / sin(angle_bot / static_cast<T>(2)));
    T dy_t = dx * (cos(angle_top / static_cast<T>(2)) / sin(angle_top / static_cast<T>(2))) + d; 

    Polygon<T> Channel({-dx, static_cast<T>(0), dx, dx, static_cast<T>(0), -dx},
                       {static_cast<T>(0), dy_b, static_cast<T>(0), d, d + dy_t, d},
                       {2, 5}, false);

    std::cout << std::setprecision(5) << std::fixed << std::showpoint << Channel << std::endl;
    return Channel;
}

// @m > 0 - Infinite horizon
// @m = 0 - Critical horizon
// @m < 0 - Finite horizon
template <typename T>
inline Polygon<T> getParallelChannel(const T &internalAngle, const T &dx, const T &m)
{
    std::cout << "Creating parallel channel with internal angle = " << toDeg(internalAngle) << " ..." << std::endl;
    T dy = dx * (cos(internalAngle / static_cast<T>(2)) / sin(internalAngle / static_cast<T>(2)));
    T d = dy + m * dy;

    if (d < 0.)
        std::cout << "Misc::getParallelChannel - Invalid pore width." << std::endl;

    Polygon<T> Channel({-dx, static_cast<T>(0), dx, dx, static_cast<T>(0), -dx},
                       {static_cast<T>(0), dy, static_cast<T>(0), d, d + dy, d},
                       {2, 5}, false);

    std::cout << std::setprecision(5) << std::fixed << std::showpoint << Channel << std::endl;
    return Channel;
}

// @m > 0 - Infinite horizon
// @m = 0 - Critical horizon
// @m < 0 - Finite horizon
template <typename T>
inline Polygon<T> getParallelChannel2(const T &dy, const T &dx, const T &d)
{
    Polygon<T> Channel({-dx, static_cast<T>(0), dx, dx, static_cast<T>(0), -dx},
                       {static_cast<T>(0), dy, static_cast<T>(0), d, d + dy, d},
                       {2, 5}, false);
    std::cout << std::setprecision(5) << std::fixed << std::showpoint << Channel << std::endl;
    return Channel;
}

template <typename T>
inline Polygon<T> getParallelChannel3(const T &alpha, const T &dx, const T &d)
{
    T dy = dx * (cos(alpha / static_cast<T>(2)) / sin(alpha / static_cast<T>(2)));
    Polygon<T> Channel({-dx, static_cast<T>(0), dx, dx, static_cast<T>(0), -dx},
                       {static_cast<T>(0), dy, static_cast<T>(0), d, d + dy, d},
                       {2, 5}, false);
    std::cout << std::setprecision(5) << std::fixed << std::showpoint << Channel << std::endl;
    return Channel;
}

template <typename T>
inline Polygon<T> getNonParallelChannel(const T &dx, const T &dyb, const T &dyt, const T &d)
{
    Polygon<T> Channel({-dx, static_cast<T>(0), dx, dx, static_cast<T>(0), -dx},
                       {static_cast<T>(0), dyb, static_cast<T>(0), d, d + dyt, d},
                       {2, 5}, false);
    std::cout << std::setprecision(5) << std::fixed << std::showpoint << Channel << std::endl;
    return Channel;
}

// @brief Generalised channel (parallel and non-parallel) determined by the 
// @isParallel true when @alpha1 == @alpha2. Use this to avoid double comparisons
// @alpha1,@alpha2 Bottom and top aperture angles
// @ratio Acts as either the pore width (explicitly) or as a ratio
template <typename T>
inline Polygon<T> getChannel(const T &dx, T &alpha1, T &alpha2, const T &ratio, const bool& isParallel)
{

    if (isParallel) 
	alpha1 = alpha2;

    T dyb = dx * cos(alpha1 / (T)2) / sin(alpha1 / (T)2);
    T dyt = dx * cos(alpha2 / (T)2) / sin(alpha2 / (T)2);
    T d = ratio * dyb;

    Polygon<T> Channel({-dx, static_cast<T>(0), dx, dx, static_cast<T>(0), -dx},
                       {static_cast<T>(0), dyb, static_cast<T>(0), d, d + dyt, d},
                       {2, 5}, false);
    
    std::cout << "Generated ";
    if (d < dyb)
	std::cout << "finite horizon ";
    else
	std::cout << "infinite horizon ";
    if (isParallel)
	std::cout << "parallel channel." << std::endl;
    else 
	std::cout << "non-parallel channel." << std::endl;
    std::cout << std::setprecision(5) << std::fixed << std::showpoint << Channel << std::endl;

    return Channel;
}


// ========================================= CHANNEL CREATORS END =========================================

// @brief Returns the indices which result in a sorted @v
// @source https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
template <typename T>
std::vector<size_t> order(const std::vector<T> &v) {

    // initialize original index locations
    std::vector<size_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values 
    std::stable_sort(idx.begin(), idx.end(),
         [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

    return idx;
}

// @brief Returns number of times a sequence of the same number is repeated
// @example v = { -1, -1, 1, 1, -1, 1, -1, -1 }
//	    return = { 2, 2, 1, 1, 2 }
//
std::vector<size_t> ConsecutiveCounts(const std::vector<int>& v) {
    size_t current_count = 1;
    std::vector<size_t> counts;
    for (size_t i = 0; i < v.size(); ++i) {
        if (v[i] == v[i+1]) {
            current_count++;
        } else {
            counts.push_back(current_count);
            current_count = 1;
        }
    }
    return counts;
}

// @brief Implementation of R's diff function with lag = differences = 1
// NOTE: @result[0] = @v[0], while @result[1] = @v[1] - @v[0]
template <typename T>
std::vector<T> diff(const std::vector<T> &v) {
    std::vector<T> result(v.size());
    std::adjacent_difference(v.begin(), v.end(), result.begin()); 
    result.erase(result.begin()); // 
    return result;

}

// @multipler - Must be an integer if T is an integer 
template <typename T>
std::vector<std::vector<T>> getMomentPowers_INT(const std::vector<int_ll>& dx, const size_t& nMoments) 
{
    std::vector<std::vector<T>> DXn(nMoments, std::vector<T>(dx.size(), -1));
    for (size_t i = 0; i < dx.size(); ++i)
    {
	T dxBase1 = (T)std::abs(dx[i]);
	T dxBase2 = dxBase1;
	for (size_t j = 0; j < nMoments; ++j)
	{
	    DXn[j][i] = dxBase1; // 
	    dxBase1 = dxBase1 * dxBase2;
	}
    }
    return DXn;
}


// @brief Returns the correlation function v[i]*v[i+m] for i = 0,1,...,v.size() and
template <typename T>
std::vector<T> Correlate(const std::vector<T>& v, const size_t& n) {
    std::vector<T> Correl(0, v.size());
    for (size_t i = 0; i < v.size() - n; ++i) 
	Correl[i] = v[i] * v[i+n];
    return Correl;
}

std::vector<int> getUniqueElements(const std::vector<int>& v) {
    std::vector<int> result = v;
    std::sort(result.begin(), result.end());
    std::vector<int>::iterator ip;
    ip = std::unique(result.begin(), result.begin() + result.size());
    result.resize(std::distance(result.begin(), ip));
    return result;
}



//// @return [vector{point}] Time (x) and step distance (y)
//template <typename T>
//std::vector<Point<T>> getLaminarSteps(const Trajectory<T>& Traj, const std::vector<int>& FlowDirs, const T& minimum_step_distance) {
//
//    std::vector<T> fps = Traj.getLiftedFreePaths();
//    std::vector<T> x = Traj.getLiftedPositions_x();
//
//    std::cout << FlowDirs << std::endl;
//    std::cout << fps << std::endl;
//    std::cout << x << std::endl;
//
//    std::vector<Point<T>> steps;
//    steps.reserve(FlowDirs.size());
//    size_t s = 1;
//    while (s < FlowDirs.size()) {
//	T current_sum_time = fps[s-1];
//	T x0 = x[s-1];
//
//	while (FlowDirs[s] == FlowDirs[s-1]) {
//	    current_sum_time += fps[s];
//	    s++;
//	    if (s == FlowDirs.size()) // Lifetime laminar stretch
//		break;
//	}
//
//	T current_sum_pos = 0;
//	if (s > 1)
//	    current_sum_pos = x[s-1] - x0;
//	else 
//	    current_sum_pos = x[0];
//	s++;
//
//	if (std::abs(current_sum_pos) > minimum_step_distance) 
//	    steps.push_back(Point<T>(current_sum_time, current_sum_pos));
//    }
//    return steps;
//}

#endif 
