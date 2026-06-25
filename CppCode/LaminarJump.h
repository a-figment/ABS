#ifndef LAMINARJUMP_H_INCLUDED
#define LAMINARJUMP_H_INCLUDED

template <typename T>
struct LaminarJump {

    LaminarJump() = default;

    // @N collisions is the maximum number of steps possible
    LaminarJump(const uint& n, const Point<T>& minstep, const T& a) : N(n), theta(a), tauxibound(minstep) {
	tau.reserve(N);
	xi.reserve(N);
	nCollisions.reserve(N); // Number of collisions along each jump
	nDirectionChanges = 1; // Number of changes in flow direction. Defaults to N
    }

    // @N collisions is the maximum number of steps possible
    LaminarJump(const uint& n, const T& a) : N(n), theta(a) {
	tauxibound = Point<T>(-1,-1);
	tau.reserve(N);
	xi.reserve(N);
	nCollisions.reserve(N); // Number of collisions along each jump
	nDirectionChanges = 1; // Number of changes in flow direction. Defaults to N
    }

    void Update_Unfiltered(const Point<T>& jump, const size_t& steps_in_flow) {
	tau.push_back(jump.x);
	xi.push_back(jump.y);
	nCollisions.push_back(steps_in_flow);
	nDirectionChanges = N;
    }

    void Update(const Point<T>& jump, const size_t& steps_in_flow) {
	// Filter out small jumps
	if (jump.x > tauxibound.x && std::abs(jump.y) > tauxibound.y) { 
	    tau.push_back(jump.x);
	    xi.push_back(jump.y);
	    nCollisions.push_back(steps_in_flow);
	    if (steps_in_flow <= N)
		nDirectionChanges++;
	} 
    }

    size_t getSize() const { return tau.size(); }
    std::vector<int> getDirectionChanges() const { 
	std::vector<int> nDC(getSize(), nDirectionChanges); 
	return nDC;
    }
    std::vector<T> getInitialVelocity() const { 
	std::vector<T> IV(getSize(), theta); 
	return IV;
    }

    std::vector<T> getJumpVelocity() const { 
	std::vector<T> jv(tau.size(), 0);
	for (size_t i = 0; i < jv.size(); ++i)
	    jv[i] = std::abs(xi[i]) / tau[i];
	return jv;
    }

    Point<T> getMinMaxJumpVelocity() const { 
	std::vector<T> jv = getJumpVelocity();
	T minjv = *min_element(jv.begin(), jv.end());
	T maxjv = *max_element(jv.begin(), jv.end());
	return Point<T>(minjv, maxjv);
    }


    int nDirectionChanges; 
    uint N;
    T theta;
    Point<T> tauxibound;
    std::vector<T> tau, xi;
    std::vector<int> nCollisions;
};

template <typename T>
std::ostream &operator<<(std::ostream &os, const LaminarJump<T> &LJ)
{
    os << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;
    std::vector<int> DC = LJ.getDirectionChanges();
    std::vector<T> IV = LJ.getInitialVelocity();
    for (size_t i = 0; i < LJ.getSize(); ++i) 
	os << LJ.tau[i] << " " << LJ.xi[i] << " " << LJ.nCollisions[i] << " " << DC[i] << " " << IV[i] << std::endl;
    return os;

}

#endif // JUMPEVENT_H_INCLUDED
