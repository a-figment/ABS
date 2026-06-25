#ifndef COARSEGRAIN_H_INCLUDED
#define COARSEGRAIN_H_INCLUDED

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <deque>
#include <numeric>
#include "Helper.h"
#include "JumpEvent.h"
#include "OnlineStatistics.h"
//#include <boost/math/statistics/linear_regression.hpp>



// @brief Rolling maximum. For rolling minimum, multiply @x by -1. 
// @x The vector to apply rolling min/max to
// @w The window half size (in indices, not true coordinates)
// @source tutorialspoint.com/sliding-window-maximum-in-cplusplus
template <typename T>
std::vector<T> rollingminmax(const std::vector<T>& x, const uint& W) {
    std::vector<T> rmm;
    std::deque<T> dq;
    uint w = 2*W + 1;

    if (w >= (uint)x.size())
	return {};

    for(uint i = 0; i < w; i++) {
	while(!dq.empty() && x[dq.back()] < x[i])
           dq.pop_back();
       dq.push_back(i);
    }
    for(uint i = w; i < x.size(); i++) {
       rmm.push_back(x[dq.front()]);
       while(!dq.empty() && dq.front() < (i - w + 1))
           dq.pop_front();
       while(!dq.empty() && x[dq.back()] < x[i])
           dq.pop_back();
       dq.push_back(i);
    }
    rmm.push_back(x[dq.front()]);
    return rmm;
}

template <typename T>
void ScaleNormalComponents(std::vector<Point<T>> CNorms, const T& tau, const T& xi) {
    for (size_t i = 0; i < CNorms.size(); ++i) {
        CNorms[i].x = CNorms[i].x / tau;
        CNorms[i].y = CNorms[i].y / std::abs(xi);
    }
}

// @brief
// @mins Result of @argminmax(-x, y, w)
// @maxs Result of @argminmax(+x, y, w)
template <typename T>
class CoarseGrain
{
    public:

	CoarseGrain(const std::vector<T>& X, const std::vector<T>& Y) {
	    x = X;
	    y = Y;
	}

	void Grain(const uint& w) {

	    std::vector<size_t> mins = argminmax(multiply<T>(y, (T)-1), w);
	    std::vector<size_t> maxs = argminmax(y, w);

	    // Signal for caller that no turning points were found
	    no_tps = false; // Assume there are tps
	    if (mins.size() == 0 || maxs.size() == 0) {
		no_tps = true;
	    } else {
	    	setTurningPointIDs(mins, maxs);
	    	setTurningPointCoords();
	    }

	}

	std::vector<size_t> argminmax(const std::vector<T> &y, const uint &w);
	void setTurningPointIDs(const std::vector<size_t>& mins, const std::vector<size_t>& maxs);
	//uint WindowSizeSelect(const uint& init_wsize, const uint& decrement, const T& R2_tol);
	uint WindowSizeSelect(const uint& init_wsize, const uint& decrement, const T& dist_thresh);

	void setTurningPointCoords() {
	    xtp.resize(idx_of_sort.size());
	    ytp.resize(idx_of_sort.size());
	    for (size_t k = 0; k < idx_of_sort.size(); ++k) {
		xtp[k] = x[idx_of_sort[k]];
		ytp[k] = y[idx_of_sort[k]];
	    }
	}

	std::vector<Point<T>> getTurningPoints() const {
	    std::vector<Point<T>> TurningPoints(idx_of_sort.size());
	    for (size_t k = 0; k < TurningPoints.size(); ++k) 
		TurningPoints[k] = Point<T>(xtp[k], ytp[k]);
	    return TurningPoints;
	}

	std::vector<T> getIDs() const { return idx_of_sort; }
	std::vector<T> getTurningPoints_x() const { return xtp; }
	std::vector<T> getTurningPoints_y() const { return ytp; }

	std::vector<T> getSuccessiveJumpTime() const { return diff<T>(xtp); }
	std::vector<T> getSuccessiveJumpDistance() const { return diff<T>(ytp); }

	size_t getNumDetections(const T& dxs) {
	    if (!isGrainable())
		return 0;
	    return getNumPruned(dxs);

	}

	// @dxs [numeric] Jump distance threshold
	// @dxt [numeric] Jump time threshold
	// @max_pd [numeric] threshold
	// @mean_pd [numeric] threshold
	JumpEvent<T> getJumpEvent(const T& dx_min, const T& dy_min, const T& mean_pd, const size_t& max_num_events) { 
	    assert(ytp.size() == xtp.size());
	    assert(xtp.size() > 1);

	    JumpEvent<T> JE(xtp.size());

	    size_t num_detections = 0;
	    for (size_t i = 0; i < xtp.size()-1; ++i) {

		T xdiff = xtp[i+1] - xtp[i]; // Jump time
		T ydiff = ytp[i+1] - ytp[i]; // Jump distance

		//T vel = ydiff / xdiff;
	    	//if ((std::abs(vel) > (std::abs(dy_min) / dx_min)) && std::abs(ydiff) > dy_min) {
	    	if (xdiff > dx_min && std::abs(ydiff) > dy_min) {

		    // Get indices between turning points
		    std::vector<size_t> indices_between((idx_of_sort[i+1] - idx_of_sort[i]) + 1);
		    std::iota(indices_between.begin(), indices_between.end(), idx_of_sort[i]);

		    // Get position of particle between turning point
		    std::vector<T> xx = index(x, indices_between);
		    std::vector<T> yy = index(y, indices_between);
		
		    Point<T> A(xtp[i], ytp[i]);
		    Point<T> B(xtp[i+1], ytp[i+1]);
		    //std::cout << A << ", " << xx[0] << " " << yy[0] << std::endl;
		    //std::cout << B << ", " << xx.back() << " " << yy.back() << std::endl;

		    std::vector<T> perp_dists = getPerpendicularDistances<T>(A, B, xx, yy);
		//  std::vector<Point<T>> norms = getPointNormals<T>(A, B, xx, yy);
	    
		//  ScaleNormalComponents(norms, xdiff, ydiff); // Divide x,y components of normal by tau and xi
		//  std::vector<T> sq(norms.size(), 0);
		//  for (size_t i = 0; i < norms.size(); ++i)
		//	sq[i] = sqrt(norms[i].x + norms[i].y);

		    //std::vector<T> mean_and_var = StatSummary<T>(perp_dists);
		    //std::vector<T> mean_and_var = StatSummary<T>(sq);

		    //T mean = mean_and_var[0];
		    //if (mean < mean_pd) {
		        //JE.Update(num_detections, xtp[i], ytp[i], xdiff, ydiff, mean, std);
		        //JE.Update(num_detections, xtp[i], ytp[i], xdiff, ydiff, mean, mean_and_var[1]);
		        JE.Update(num_detections, xtp[i], ytp[i], xdiff, ydiff, 0, 0);
		        num_detections++;

			if (num_detections > max_num_events)
			    return JE;
		    //}
	    	}
	    }
	    return JE;
	}


	uint getNumPruned(const T& dxs, const T& dxt) {
	    std::vector<T> jd = getSuccessiveJumpDistance();
	    std::vector<T> jt = getSuccessiveJumpTime();
	    uint count = 0;
	    for (size_t i = 0; i < jd.size(); ++i) 
		if (std::abs(jd[i]) > dxs && jt[i] > dxt) 
		    count++;
	    return count;
	}

	std::vector<T> getPairwiseJumpTime() const { 
	    std::vector<T> jt;
	    jt.reserve(xtp.size());
	    for (size_t i = 0; i < xtp.size()-1; i+=2) 
		jt.push_back(std::abs(xtp[i] - xtp[i+1]));
	    return jt;
	}

	std::vector<T> getPairwiseJumpDistance() const {
	    std::vector<T> jd;
	    jd.reserve(xtp.size());
	    for (size_t i = 0; i < ytp.size()-1; i+=2) 
		jd.push_back(ytp[i] - ytp[i+1]);
	    return jd;
	}

	bool isGrainable() const { return (!no_tps); } // Grainable when there are tp's

    private:
	std::vector<T> x;
	std::vector<T> y;
	//uint w;

	bool no_tps;

	std::vector<size_t> idx_of_sort; // Defines indices of critical points in correct order. 
					 // Used to index @x, @y. e.g. jump_time = diff(x[idx_of_sort])
	// x-y coordinate of turning points
	std::vector<T> xtp;
	std::vector<T> ytp;


};

// @brief Central function of the class. Re-maps indices of turning points found by @argminmax into indices which locate in true coordinates the beginning/end of a jump
template <typename T>
void CoarseGrain<T>::setTurningPointIDs(const std::vector<size_t>& mins, const std::vector<size_t>& maxs) {

    // x-coord of tp's
    std::vector<T> t_mins(mins.size());
    std::vector<T> t_maxs(maxs.size());
    for (size_t i = 0; i < t_mins.size(); ++i) 
        t_mins[i] = x[mins[i]];
    for (size_t i = 0; i < t_maxs.size(); ++i) 
        t_maxs[i] = x[maxs[i]];
    
    // Combine min/max tp's and get the indices of sorted values
    std::vector<T> taus = t_mins;
    taus.insert(taus.end(), t_maxs.begin(), t_maxs.end());
    std::vector<size_t> tau_ids = order(taus);
    
    // Combine and sort turning point locations
    std::vector<size_t> min_max_idx = mins;
    min_max_idx.insert(min_max_idx.end(), maxs.begin(), maxs.end());
    
    idx_of_sort.resize(tau_ids.size());
    for (size_t i = 0; i < tau_ids.size(); ++i) 
	idx_of_sort[i] = min_max_idx[tau_ids[i]];
}

// @brief 
template <typename T>
std::vector<size_t> CoarseGrain<T>::argminmax(const std::vector<T> &y, const uint &w)
{
    std::vector<T> minmax = rollingminmax<T>(y, w);
    if (minmax.size() == 0)
	return {};

    std::vector<size_t> ids;
    for (size_t i = w; i < y.size()-w; ++i) {
	T delta = minmax[i-w] - y[i];
	if (delta <= 0) {
	    ids.push_back(i);
	}
    }
    return ids;
}




// @brief The straight lines passing through the turning points  AB = (@xtp, @ytp) 
//	  examined for the 
// @D is the acceptance tolerance for maximum deviation from the line AB
// @decrement is the amount to decrement from @w upon a failed acceptance test
// @return > 0 The desired window size. 
// @return = 0 No size was found to satisfy all jumps
template <typename T>
uint CoarseGrain<T>::WindowSizeSelect(const uint& init_wsize, const uint& decrement, const T& dist_thresh) {

    //using boost::math::statistics::simple_ordinary_least_squares;

    size_t prev_idx_size = 1;
    uint ideal_size = init_wsize;
    for (uint w = init_wsize; w > decrement; w -= decrement) {
	//std::cout << "Window size = " << w << std::endl;
	
	Grain(w); // Sets the turning points
	//std::cout << "Num. T.P. = " << idx_of_sort.size() << std::endl;

	size_t num_ideal = 0; // 
	if (isGrainable()) {
	    if(2*num_ideal != prev_idx_size) { // Re-try for next window size
		prev_idx_size = idx_of_sort.size();

    		for (size_t i = 0; i < idx_of_sort.size()-1; ++i) { // Consider the lines between all turning points
    	
		  // Start and end points of a detected jump
    		  T P1x = x[idx_of_sort[i]];
    		  T P1y = y[idx_of_sort[i]];
    		  T P2x = x[idx_of_sort[i + 1]];
    		  T P2y = y[idx_of_sort[i + 1]];
    		  T den = sqrt((P2x - P1x)*(P2x - P1x) + (P2y - P1y)*(P2y - P1y));
    	
    		  size_t sz = idx_of_sort[i+1] - idx_of_sort[i];
    		  std::vector<T> perp_dists(sz,0); // 
    		  for (size_t j = idx_of_sort[i]; j < idx_of_sort[i+1]; ++j) { // Construct the line (X,Y) between current TP's
    		      // Raw data points between turning point coords
    		      T x0 = x[j];
    		      T y0 = y[j];
    	
    		      T num = std::abs((P2x - P1x) * (P1y - y0) - (P1x - x0)*(P2y - P1y));
    		      perp_dists[j - idx_of_sort[i]] = num / den;
    		  }
    	

    		  //T mean = std::reduce(perp_dists.begin(), perp_dists.end()) / (T)perp_dists.size();
		  //T mean = std::accumulate( v.begin(), v.end(), 0.0)/v.size();        
    		  T M = *std::max_element(perp_dists.begin(), perp_dists.end());

		  //std::cout << "Mean = " << mean << ", max = " << M << std::endl;
    		  
    		  if (M > dist_thresh) { // The line between this pair of TP's failed test. Re-try with a new size
    		      break;
    		  } else {
		      std::cout << w << std::endl;
    		      ideal_size = w;
    		      num_ideal++;
    		  }
		    
		}
	    } else {
		break;
	    }
	}
    }

    if (!isGrainable() || ideal_size == init_wsize)
	return 0;
    return ideal_size;

}



#endif 
