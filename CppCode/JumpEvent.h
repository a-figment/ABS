#ifndef JUMPEVENT_H_INCLUDED
#define JUMPEVENT_H_INCLUDED

template <typename T>
struct JumpEvent {

    JumpEvent() = default;
    JumpEvent(const size_t& res) {
    	count.reserve(res);
	t0.reserve(res);
    	x0.reserve(res);
    	tau.reserve(res);
    	xi.reserve(res);
    	mean_pd.reserve(res);
    	var_pd.reserve(res);
    	//numLR.reserve(res);
    }

    void Update(const size_t& cnt, 
	   const T& t_st,
	   const T& x_st,
	   const T& jt,
	   const T& jd,
	   const T& mpd,
	   const T& vpd) {
	   //const Point<int>& numlr) {

    	count.push_back(cnt);
	t0.push_back(t_st);
    	x0.push_back(x_st);
    	tau.push_back(jt);
    	xi.push_back(jd);
    	mean_pd.push_back(mpd);
    	var_pd.push_back(vpd);
    	//numLR.push_back(numlr);
    }

    void setFirstKJumps(const size_t& KK) {
	size_t res = KK;

	std::vector<size_t> _count;
	std::vector<T> _t0, _x0, _tau, _xi, _mean_pd, _var_pd;

	_count.reserve(res);
	_t0.reserve(res);
    	_x0.reserve(res);
    	_tau.reserve(res);
    	_xi.reserve(res);
    	_mean_pd.reserve(res);
    	_var_pd.reserve(res);

	int i = 0;
	while (i <= KK && i < xi.size()) {
	   _count.push_back(count[i]);
	   _t0.push_back(t0[i]);
    	   _x0.push_back(x0[i]);
    	   _tau.push_back(tau[i]);
    	   _xi.push_back(xi[i]);
    	   _mean_pd.push_back(mean_pd[i]);
    	   _var_pd.push_back(var_pd[i]);
	   i++;
	}

	count = _count;
	t0 = _t0;
	x0 = _x0;
	tau = _tau;
	xi = _xi;
	mean_pd = _mean_pd;
	var_pd = _var_pd;
    }

    // MUST be run after @setFirstKJumps
    void VelocityFilter(const T& vel, const T& eps) {
	size_t res = tau.size(); 

	std::vector<size_t> _count;
	std::vector<T> _t0, _x0, _tau, _xi, _mean_pd, _var_pd;

	_count.reserve(res);
	_t0.reserve(res);
    	_x0.reserve(res);
    	_tau.reserve(res);
    	_xi.reserve(res);
    	_mean_pd.reserve(res);
    	_var_pd.reserve(res);

	for (size_t i = 0; i < xi.size(); ++i) {
	    T v = std::abs(xi[i] / tau[i]);
	    //if ((v < vel - eps) || (v > vel + eps)) {
	    if ((v > vel - eps) && (v < vel + eps)) {
		_count.push_back(count[i]);
		_t0.push_back(t0[i]);
		_x0.push_back(x0[i]);
		_tau.push_back(tau[i]);
		_xi.push_back(xi[i]);
		_mean_pd.push_back(mean_pd[i]);
		_var_pd.push_back(var_pd[i]);
	    }
	}

	count = _count;
	t0 = _t0;
	x0 = _x0;
	tau = _tau;
	xi = _xi;
	mean_pd = _mean_pd;
	var_pd = _var_pd;
    }



    std::vector<size_t> count;
    std::vector<T> t0, x0, tau, xi, mean_pd, var_pd;
    //std::vector<Point<int>> numLR;
};

template <typename T>
std::ostream &operator<<(std::ostream &os, const JumpEvent<T> &JE)
{
    os << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;
    for (size_t i = 0; i < JE.count.size(); ++i) {
	os << JE.count[i] << " " << JE.t0[i] << " " << JE.x0[i] << " " << JE.tau[i] << " " << JE.xi[i] << " " << JE.mean_pd[i] << " " << JE.var_pd[i] << std::endl;
	//" " << JE.numLR[i].x << " " << JE.numLR[i].y << std::endl;
    }
    return os;

}

#endif // JUMPEVENT_H_INCLUDED
