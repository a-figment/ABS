#ifndef SCATTER_H_INCLUDED
#define SCATTER_H_INCLUDED

const float_ldub PI = boost::math::constants::pi<float_ldub>();
const float_ldub P2 = PI/(float_ldub)2; 
const float_ldub A0 = 3*P2;
const float_ldub A20 = P2;
const float_ldub dx2 = (float_ldub)1; // 2*\delta x = 1
const float_ldub oh = 1/(float_ldub)2; // 1/2
const float_ldub dx = oh;
const float_ldub dy = oh;

// Variables of the scattering map
template <typename T>
struct ScatterVar {

    ScatterVar(const T& height, const T& angle, const T& tau, const int_ll& position, const int& bid) {
	Set(height, angle, tau, position, bid);
    }

    void Set(const T& height, const T& angle, const T& tau, const int_ll& position, const int& bid) {
	Height = height;
	Angle = angle;
	Tau = tau;
	Position = position;
	BranchID = bid;
    }

    void Print() {
	std::cout << "Angle           = " << Angle<< " (" << Angle * 180/PI << " degrees)" << std::endl;
	std::cout << "Height          = " << Height << std::endl;
	std::cout << "Position        = " << Position << std::endl;
	std::cout << "Time            = " << Tau << std::endl;
	std::cout << "Branch ID       = " << BranchID << std::endl;
    }

    T Height;
    T Angle;
    T Tau;
    int_ll Position;
    int BranchID;
};


template <typename T>
struct ScatterTraj {

    ScatterTraj() = default;

    ScatterTraj(const size_t& N) {
	idx = 0;
	Height.resize(N);
	Angle.resize(N);
	Tau.resize(N);
	Position.resize(N);
	BranchID.resize(N);
	//Update(height, angle, 0, position, -1);
    }

    void Update(const T& height, const T& angle, const T& tau, const int_ll& position, const int& bid) {
	Height[idx] = height;
	Angle[idx] = angle;
	Tau[idx] = tau;
	Position[idx] = position;
	BranchID[idx] = bid;
	idx++;
    }

    void Update(const ScatterVar<T>& SV) {
	Update(SV.Height, SV.Angle, SV.Tau, SV.Position, SV.BranchID);
    }

    std::vector<int_ll> getPositions() const { return Position; }
    std::vector<int> getBranchIDs() const { return BranchID; }
    std::vector<T> getHeights() const { return Height; }
    std::vector<T> getAngles() const { return Angle; }
    std::vector<T> getTaus() const { return Tau; }

    std::vector<int> BranchID;		// A number between 0 - 19 indicating which 
    std::vector<int_ll> Position;	// Lifted position of particle
    std::vector<T> Height, Angle, Tau;  // Entrance heights, angles and jump lengths (cumulative)
    size_t idx;
};


// # Gamma 0 MONOTONIC ANGLES
//template <typename T> 
//T a0(const T& d, const T& h0)  { return (-PI/(T)2 + 2*PI); }
template <typename T> 
T a1(const T& d, const T& h0)  { return (-atan(h0/d)); }
template <typename T> 
// We define this so that it is valid for IH, but there is no such andlge between a2 and a3
// when d = 1. 
T a2(const T& d, const T& h0)  { return (-atan((h0 + d - oh)/(dx2 + d - oh))); }
template <typename T> 
T a3(const T& d, const T& h0)  { return (-atan(h0/(dx2 + d))); }
template <typename T> 
T a4(const T& d, const T& h0)  { return (-atan(h0/(dx2 + 2*d))); }
template <typename T> 
T a5(const T& d, const T& h0)  { return atan((d-h0)/(dx2 + 2*d)); }
template <typename T> 
T a6(const T& d, const T& h0)  { return atan((dy-h0)/(dx2 + 2*d - oh)); }
template <typename T> 
T a7(const T& d, const T& h0)  { return atan((dy-h0)/(dx2 + d - oh)); }
template <typename T> 
T a8(const T& d, const T& h0)  { return atan((d-h0)/dy); } // Angle to bottom middle vertex (dx, dy)

template <typename T> 
T phi0(const T& d, const T& h0)  { return (-atan((h0+d)/(4*d))); }
template <typename T> 
T phi1(const T& d, const T& h0)  { return (-atan((h0+d)/(5*d))); }
template <typename T> 
T phi2(const T& d, const T& h0)  { return (atan((dx2 - h0)/(dx2 + 2*d))); }
template <typename T> 
T phi3(const T& d, const T& h0)  { return (atan((dx2 - h0)/(dx2 + d))); }

// # Gamma 2 angles
template <typename T> 
T a12(const T& d, const T& h0) { return atan((d - h0)/d); } // When d = 1/2 this is the angle to the center vertex otherise it is the 
						// angle to the top right vertex.
template <typename T> 
T a13(const T& d, const T& h0) { return atan((d + dy - h0)/(dx2 + d - oh)); }
template <typename T> 
T a14(const T& d, const T& h0) { return atan((d + dy + dx - h0)/(dx2 + d)); }
template <typename T> 
T a15(const T& d, const T& h0) { return atan((d + dx2 - h0)/dx2); }
template <typename T> 
T a16(const T& d, const T& h0) { return atan((d + dy - h0) / dx); }

// # Gamma 3 angles
template <typename T> 
T a17(const T& d, const T& h0) { return atan((d + dy - h0) / dx); }
template <typename T> 
T a18(const T& d, const T& h0) { return atan((2*d + dx2 - h0)/(dx2)); }
template <typename T> 
T a19(const T& d, const T& h0) { return atan((2*dy + 2*d - dx - h0)/dx); }
//template <typename T> 
//T a20(const T& d, const T& h0) { return PI/(T)2; }

template <typename T> 
int_ll kappa(const T& d, const T& h0, const T& t0) { 
    return ceil((d - h0 - tan(t0)) / (d*(tan(t0) - (T)1)));
}

template <typename T> 
int_ll gamma(const T& d, const T& h0, const T& t0) { 
    return ceil((h0 - (T)1 - d) / (d*((T)1 - tan(t0))));
}

template <typename T> 
int_ll alpha(const T& d, const T& h0, const T& t0) { 
    return ceil((-h0 - tan(t0)) / (d * ((T)1 + tan(t0))));
}

template <typename T>
class Scatterer {

    public:
	Scatterer(const T& D, const T& H0, const T& T0, const int_ll& X0, const int& BID) {

	    // Initial conditions
	    d = D;
	    nIterations = 0;
	    Assign(H0, T0, X0, 0, BID);
	}

	// Updates the field variables
	void Assign(const T& height, const T& theta, const int_ll& pos, const T& time, const int& branch) {
	    H1 = height;
	    T1 = theta;
	    X1 = pos;
	    TAU = time;
	    BranchID = branch;
	}
	
	// @brief Updates the position of a particle based on its exit angle 
	void UpdatePosition(const T& T1, int_ll& x0, const bool& sym) {
	    if (sym && (cos(T1) >= 0)) { 
		x0 = x0 - 1;
	    } else if (!sym && (cos(T1) > 0)) {
		x0 = x0 + 1;
	    } else {
		x0 = x0;
	    }
	}

	void Update(T& h0, T& t0, int_ll& x0, T& tau, int& bid, bool& sym)  {

	    if (t0 > atan((oh-h0)/oh) && t0 < atan(d - h0)) {		    // RIGHT EXIT
		Er(h0, t0, tau, bid);
		UpdatePosition(t0, x0, sym);
	    } else if ((t0 > -P2) && t0 < atan((dy-h0)/dx)) {		    // BOTTOM LEFT
		Sg0(h0, t0, tau, bid);
		UpdatePosition(t0, x0, sym);
	    } else if (t0 > atan((d + dy - h0) / dx) && t0 < P2) {	    // TOP LEFT
		Sg3(h0, t0, tau, bid);
		UpdatePosition(t0, x0, sym);
	//    } else if (d > 1 && t0 > atan((h - oh)/oh) && t0 < atan(h)) {   // BOTTOM RIGHT
	//	Sg1(h0, t0, tau, bid);
	//	UpdatePosition(t0, x0, sym);
	    } else {							    // TOP RIGHT
		Sg2(h0, t0, tau, bid);
		UpdatePosition(t0, x0, sym);
	    }
	    Assign(h0, t0, x0, tau, bid);
	}


	// @brief Iterates the scattering map one time WITHOUT RECURSION
	void Evolve_(T& h0, T& t0, int_ll& x0, T& tau, int& bid, bool& sym) {

	    if (h0 > d || h0 < 0) {
		std::cout << "h0 is out of bounds." << std::endl;
	    }

	    if ((t0 > (-P2)) && (t0 < P2)) { // Rightward moving
		Update(h0, t0, x0, tau, bid, sym);
		nIterations++;
	    } else {
		sym = true;	
		T temp = PI - t0; // Reference variable
		Update(h0, temp, x0, tau, bid, sym);
		t0 = PI - temp; // Revert 
	    }
	}

	///void Sg1(T& h0, T& t0, T& tau, int& bid) {
	///    tau = //1/cos(t0);
	///    h0  = //tan(t0) + h0;
	///    t0  = //t0;
	///    bid = //20;
	///}

	// For infinite horizon, we can exit 
	void Er(T& h0, T& t0, T& tau, int& bid) {
	    tau = 1/cos(t0);
	    h0 = tan(t0) + h0;
	    t0 = t0;
	    bid = 19;
	}

	// @brief Returns index of the Sg0 branch
	void Sg0(T& h0, T& t0, T& tau, int& bid) {
	    if (t0 > -P2 && t0 < a1(d,h0)) {
		tau = - h0 / sin(t0);
		h0 = -h0 * cot(t0);
		t0 = P2 - t0;
		bid = 0;
	    } else if (t0 > a1(d,h0) && t0 < a2(d,h0)) {
		if (t0 < phi0(d,h0)) {
		    int_ll ALPHA = alpha(d,h0,t0); // Save a calculation
		    tau = d*(alpha(d,h0,t0)+3) / cos(t0);
		    h0 = h0 + d*(ALPHA + 3)*tan(t0) + d*(ALPHA + 1);
		    t0 = t0;
		    bid = 1;
		} else if (t0 > phi0(d,h0) && t0 < phi1(d,h0)) {
		    tau = -(h0 + d) / sin(t0);
		    h0 = (h0+d)*cot(t0) + 5*d;
		    t0 = -(t0 + P2);
		    bid = 2;
		} else { // t0 > phi1
		    tau = (3*dx + 2*d)/cos(t0);
		    h0 = h0 + d*(1 + 5*tan(t0));
		    t0 = t0;
		    bid = 3;
		}
	    } else if (t0 > a2(d,h0) && t0 < a3(d,h0)) {
		tau = (d + dx2)/cos(t0);
		//h0 = h0+(1+d)*tan(t0) + d;
		h0 = h0 + (d + dx2) * tan(t0) + d;
		t0 = t0;
		bid = 4;
	    } else if (t0 > a3(d,h0) && t0 < a4(d,h0)) {
		tau = -h0/sin(t0);
		//h0 = h0*cot(t0) + 4*d;
		h0 = h0*cot(t0) + 2*d + dx2;
		t0 = -(t0 + P2);
		bid = 5;
	    } else if ((d < 1 && (t0 > a4(d,h0) && t0 < a5(d,h0)))
		    || (d > oh && (t0 > a4(d,h0) && t0 < a6(d,h0)))) {
		//tau = 4*d/cos(t0);
		tau = (dx2 + 2*d)/cos(t0);
		h0 = h0 + (dx2 + 2*d)*tan(t0);
		t0 = t0;
		bid = 6;
	    } else if (t0 > a5(d,h0) && t0 < a6(d,h0)) {
		//tau = 5*d/cos(t0);
		tau = (3*dx + 2*d)/cos(t0);
		h0 = h0 + (2 + d)*tan(t0) - d;
		t0 = t0;
		bid = 7;
	    } else if (t0 > a6(d,h0) && t0 < a7(d,h0)) {
		//tau = 4*d/cos(t0);
		tau = (dx2 + 2*d)/cos(t0);
		//h0 = -h0 + 1 - 2*tan(t0);
		h0 = dx2 - (h0 + (dx2 + 2*d) * tan(t0));
		t0 = t0 + PI;
		bid = 8;
	    } else if (t0 > a7(d,h0)  && t0 < a8(d,h0)) {
		if (t0 < phi2(d,h0)) {
		    tau = (dx2 + 2*d)/cos(t0);
		    h0 = dx2 - (h0 + (dx2 + 2*d)*tan(t0));
		    t0 = t0 + PI;
		    bid = 9;
		} else if (t0 > phi2(d,h0) && t0 < phi3(d,h0)) {
		    tau = (dx2-h0)/sin(t0);
		    h0 = dx2 + 2*d - (dx2 - h0)/tan(t0);
		    t0 = 3*P2 - t0;
		    bid = 10;
		} else { // t0 > phi3
		    tau = (dx2 + d)/cos(t0);
		    h0 = d + dx2 - (h0 + (dx2 + d)*tan(t0));
		    t0 = t0 + PI;
		    bid = 11;
		}
	    } else {
		std::cout << "Impossble angle in Sg0" << std::endl;
		std::cout << "h0 = " << h0 << std::endl;
		std::cout << "t0 = " << t0*180/PI << std::endl;
		tau = -1;
		h0 = -1;
		t0 = 0;
		bid = -1;
	    }
	}

	// @brief Returns index of the Sg2 branch
	void Sg2(T& h0, T& t0, T& tau, int& bid) {
	    if (t0 > a12(d,h0) && t0 < a13(d,h0)) {
		tau = (dx2+kappa(d,h0,t0)*d)/cos(t0);
		//h0 = h0 + d*kappa(d,h0,t0) * (tan(t0) - 1) + tan(t0);
		h0 = (h0 + (dx2 + kappa(d,h0,t0)*d)*tan(t0)) - d*kappa(d,h0,t0);
		t0 = t0;
		bid = 12;
	    } else if (t0 > a13(d,h0) && t0 < a14(d,h0)) {
		//tau = 3*d/cos(t0);
		//h0 = -h0 + (1 - tan(t0))*(1 + d);
		tau = (dx2 + d)/cos(t0);
		h0 = d + dx + dy - (h0 + (dx2 + d)*tan(t0));
		t0 = t0 + PI;
		bid = 13;
	    } else if (t0 > a14(d,h0) && t0 < a15(d,h0)) {
		//tau = (3*d - h0)/sin(t0);
		//h0 = cot(t0) * (h0 + (d + 1)*(tan(t0) - 1));
		tau = (d + dx2 - h0)/sin(t0);
		h0 = (dx2 + d) - cot(t0) * (d + dx2 - h0);
		t0 = -t0 + 3*P2;
		bid = 14;
	    } else if (t0 > a15(d,h0) && t0 < a16(d,h0)) {
		tau = dx2/cos(t0);
		h0 = 2*dy + 2*d - (h0 + tan(t0));
		t0 = t0 + PI;
		bid = 15;
	    } else {
		std::cout << "Impossble angle in Sg2" << std::endl;
		std::cout << "h0 = " << h0 << std::endl;
		std::cout << "t0 = " << t0*180/PI << std::endl;
		tau = -1;
		h0 = -1;
		t0 = 0;
		bid = -1;
	    }
	}

	// @brief Returns index of the Sg3 branch
	void Sg3(T& h0, T& t0, T& tau, int& bid) {
	    if (t0 > a17(d,h0) && t0 < a18(d,h0)) {
		tau = dx2/cos(t0);
		h0 = 2*dy + 2*d - (h0 + tan(t0));
		t0 = t0 + PI;
		bid = 16;
	    } else if (t0 > a18(d,h0) && t0 < a19(d,h0)) {
		tau = (2*dy + 2*d - h0)/sin(t0);
		h0 = dx2 - cot(t0)*(2*dy + 2*d - h0);
		t0 = (3*P2) - t0;
		bid = 17;
	    } else if (t0 > a19(d,h0) && t0 < P2) { // Upper bound is a20
		tau = (dx2 + d*(1 + gamma(d,h0,t0)) - h0)/sin(t0);
		//h0 = cot(t0) * (-h0 + 1 + d*gamma(d,h0,t0)*(1-tan(t0)) + d*(1 + tan(t0)));
		h0 = cot(t0) * (dx2 + d*(1 + gamma(d,h0,t0)) - h0) - d*(gamma(d,h0,t0) - 1);
		t0 = P2 - t0;
		bid = 18;
	    } else {
		std::cout << "Impossble angle in Sg3" << std::endl;
		std::cout << "h0 = " << h0 << std::endl;
		std::cout << "t0 = " << t0*180/PI << std::endl;
		tau = -1;
		h0 = -1;
		t0 = 0;
		bid = -1;
	    }
	}

	T getHeight() const { return H1; }
	T getAngle() const { return T1; }
	int_ll getPosition() const { return X1; }
	T getTime() const { return TAU; }
	int getBranchID() const { return BranchID; }
	size_t getIterations() const { return nIterations; }

    private:
	T d;
	size_t nIterations;
	int BranchID;    // A number between 0 - 19 indicating which 
	T H1, T1, TAU;
	int_ll X1;

};

template <typename T>
std::ostream &operator<<(std::ostream &os, const Scatterer<T> &S)
{
  os << "Angle           = " << S.getAngle() << " (" << S.getAngle() * 180/PI << " degrees)" << std::endl;
  os << "Height          = " << S.getHeight() << std::endl;
  os << "Position        = " << S.getPosition() << std::endl;
  os << "Time            = " << S.getTime() << std::endl;
  os << "Branch ID       = " << S.getBranchID() << std::endl;
  os << "Num. Iterations = " << S.getIterations() << std::endl;
  return os;
}

// @ID The branch ID from which to sample a point
// @TODO: What happens if we collect also the crossing time corresponding to each interval?
// TODO: The problem here is that we can only check angles on the right half plane. It needs to be extended to 
//	 return points also on the left half plane giving angles in [0, 2pi). Further, we need to find the area
//	 of each region. Otherwise, this method is much faster than the rejection sampling approach
//
template <typename T>
Point<T> getPointInBranch_(const int& ID, const T& eps, const T& d) {

    Random<T> URA, URH;
    T A = -1;
    T H = -1;

    int bid = -1;
    while (bid != ID) {
	H = URH.getUniformRandom(eps, d-eps);
	A = URA.getUniformRandom(-P2, P2);

	if (A > -P2 && A < atan((dy - H) / dx)) {  // BOTTOM LEFT
	    if (A > -P2 && A < a1(d,H)) {
		bid = 0;
	    } else if (A > a1(d,H) && A < a2(d,H)) {
		if (A < phi0(d,H)) {
		    bid = 1;
		} else if (A > phi0(d,H) && A < phi1(d,H)) {
		    bid = 2;
		} else { // A > phi1
		    bid = 3;
		}
	    } else if (A > a2(d,H) && A < a3(d,H)) {
		bid = 4;
	    } else if (A > a3(d,H) && A < a4(d,H)) {
		bid = 5;
	    } else if ((d < 1 && (A > a4(d,H) && A < a5(d,H)))
		    || (d > oh && (A > a4(d,H) && A < a6(d,H)))) {
		    bid = 6;
	    } else if (A > a5(d,H) && A < a6(d,H)) {
		bid = 7;
	    } else if (A > a6(d,H) && A < a7(d,H)) {
		bid = 8;
	    } else if (A > a7(d,H) && A < a8(d,H)) {
		if (A < phi2(d,H)) {
		    bid = 9;
		} else if (A > phi2(d,H) && A < phi3(d,H)) {
		    bid = 10;
		} else { // A > phi3
		    bid = 11;
		}
	    } else {
		std::cout << "Angle in THETA_g0" << std::endl;
		A = -1;
		H = -1;
	    }
	}

	//if (((A > atan((dy - H) / dx) && d < 1) || (A > atan((d - H)) && d > 0.5)) && A < atan((dy + d - H)/dx)) {  // TOP RIGHT
	if (A > a12(d,H) && A < atan((d + dy - H)/dx)) {
	    if (A > a12(d,H) && A < a13(d,H)) {
		bid = 12;
	    } else if (A > a13(d,H) && A < a14(d,H)) {
		bid = 13;
	    } else if (A > a14(d,H) && A < a15(d,H)) {
		bid = 14;
	    } else if (A > a15(d,H) && A < a16(d,H)) {
		bid = 15;
	    } else {
		std::cout << "Angle in THETA_g2" << std::endl;
		A = -1;
		H = -1;
	    }
	}

	if (A > atan((dy + d - H)/dx) && A < P2) { // TOP LEFT
	    if (A > a17(d,H) && A < a18(d,H)) {  
		bid = 16;
	    } else if (A > a18(d,H) && A < a19(d,H)) {
		bid = 17;
	    } else if (A > a19(d,H) && A < P2) { // Upper bound is a20
		bid = 18;
	    } else {
		std::cout << "Angle in THETA_g3" << std::endl;
		A = -1;
		H = -1;
	    }
	}

	if (d > 0.5 && (A > atan((dy - H)/dx) && A < atan((d - H)/dx2))) { // EXIT RIGHT Er
	    bid = 19;
	}
    }
    return Point<T>(H,A);
}


template <typename T>
// @ID The branch ID from which to sample a point
std::pair<std::vector<T>, std::vector<T>> getPointsInBranch_(const int& ID, const T& eps, const T& d, const size_t& nPoints) {

    std::vector<T> Heights(nPoints, 0);
    std::vector<T> Angles(nPoints, 0);
    for (size_t i = 0; i < nPoints; ++i) {
	Point<T> P = getPointInBranch_(ID, eps, d);
	Heights[i] = P.x;
	Angles[i] = P.y;
    }
    return std::make_pair(Heights, Angles);
}

// @ID Check if the point (h, t) is in the region @ID 
template <typename T>
bool isPointInBranch(const T& H, const T& a, const std::vector<int>& ID, const T& d) {

    T A = a;
    if (A > P2)
	A = PI - a;

    int bid = -1;
    if (A > -P2 && A < atan((dy - H) / dx)) {  // BOTTOM LEFT
	if (A > -P2 && A < a1(d,H)) {
	    bid = 0;
	} else if (A > a1(d,H) && A < a2(d,H)) {
	    if (A < phi0(d,H)) {
		bid = 1;
	    } else if (A > phi0(d,H) && A < phi1(d,H)) {
		bid = 2;
	    } else { // A > phi1
		bid = 3;
	    }
	} else if (A > a2(d,H) && A < a3(d,H)) {
	    bid = 4;
	} else if (A > a3(d,H) && A < a4(d,H)) {
	    bid = 5;
	} else if ((d < 1 && (A > a4(d,H) && A < a5(d,H)))
		|| (d > oh && (A > a4(d,H) && A < a6(d,H)))) {
	    bid = 6;
	} else if (A > a5(d,H) && A < a6(d,H)) {
	    bid = 7;
	} else if (A > a6(d,H) && A < a7(d,H)) {
	    bid = 8;
	} else if (A > a7(d,H) && A < a8(d,H)) {
	    if (A < phi2(d,H)) {
		bid = 9;
	    } else if (A > phi2(d,H) && A < phi3(d,H)) {
		bid = 10;
	    } else { // A > phi3
		bid = 11;
	    }
	} else {
	    std::cout << "Angle in THETA_g0" << std::endl;
	}
    }

    if (A > a12(d,H) && A < atan((d + dy - H)/dx)) {
	if (A > a12(d,H) && A < a13(d,H)) {
	    bid = 12;
	} else if (A > a13(d,H) && A < a14(d,H)) {
	    bid = 13;
	} else if (A > a14(d,H) && A < a15(d,H)) {
	    bid = 14;
	} else if (A > a15(d,H) && A < a16(d,H)) {
	    bid = 15;
	} else {
	    std::cout << "Angle in THETA_g2" << std::endl;
	}
    }

    if (A > atan((dy + d - H)/dx) && A < P2) { // TOP LEFT
	if (A > a17(d,H) && A < a18(d,H)) {  
	    bid = 16;
	} else if (A > a18(d,H) && A < a19(d,H)) {
	    bid = 17;
	} else if (A > a19(d,H) && A < P2) { // Upper bound is a20
	    bid = 18;
	} else {
	    std::cout << "Angle in THETA_g3" << std::endl;
	}
    }

    if (d > 0.5 && (A > atan((dy - H)/dx) && A < atan((d - H)/dx2))) { // EXIT RIGHT Er
	bid = 19;
    }
    return (In<int>(bid,ID));
}


// Rejection samples uniformly distributed points until a point in branch @ID is found
template <typename T>
std::pair<std::vector<T>, std::vector<T>> getPointsInBranch_Rejection(const std::vector<int>& ID, const T& eps, const T& d, const size_t& nPoints) {

    std::vector<T> Heights(nPoints, 0);
    std::vector<T> Angles(nPoints, 0);
    Random<T> URHeights, URAngles;
    for (size_t i = 0; i < nPoints; ++i) {

	T hTest = URHeights.getUniformRandom(eps, d - eps);
	T aTest = URAngles.getUniformRandom(-P2+eps, 3*P2 - eps);

	while (!isPointInBranch(hTest, aTest, ID, d)) {
	    hTest = URHeights.getUniformRandom(eps, d - eps);
	    aTest = URAngles.getUniformRandom(-P2+eps, 3*P2 - eps);
	}

	Heights[i] = hTest;
	Angles[i] = aTest;
    }
    return std::make_pair(Heights, Angles);
}

template <typename T>
ScatterTraj<T> getTrajectory(const T& Height, const T& Angle, const T& d, const std::vector<T> Time) {

    T dt = Time[3] - Time[2]; // Bit dodgy
    // IC's
    T h0 = Height;
    T T0 = Angle;
    T tau = 0;
    int_ll x0 = 0;
    int bid = -1;
    bool sym = false;

    // Evolve in continuous time
    size_t tidx = 0;
    T TAU = tau;  // Current time
    T pTAU = TAU; // Previous time

    Scatterer<T> S(d, h0, T0, x0, bid);	// Evolves the map
    ScatterVar<T> SV(h0, T0, pTAU, x0, bid);      // Variables of the map
    ScatterTraj<T> ST(Time.size()); // Stores variables 

    // Evolve map until time exceeds target time
    while (TAU < Time.back()) {

	S.Evolve_(h0, T0, x0, tau, bid, sym);
	sym = false;
	TAU += tau;

	if (tau > 6*dt) {
	    std::cout << " ===================================== " << std::endl;
	    std::cout << "Crossing time tau = " << tau << 
	    " was longer than step time = CHECK THIS!!!" << 6*dt << std::endl;
	    std::cout << S << std::endl;
	    SV.Print();
	    std::cout << " ===================================== " << std::endl << std::endl;
	    break;
	}

	if (tidx < Time.size() && TAU >= Time[tidx]) { // Exceeded a target time
	    // Is current iterate closer to Time than previous iterate?
	    if (abs(TAU - Time[tidx]) < abs(pTAU - Time[tidx])) {
		ST.Update(h0, T0, TAU, x0, bid);
	    } else {
		ST.Update(SV);
	    }
	    tidx++;
	}
	pTAU = TAU;
	SV.Set(h0, T0, pTAU, x0, bid); // Previous state
    }

    return ST;
}

template <typename T>
ScatterTraj<T> getTrajectory(const T& Height, const T& Angle, const T& d, const std::vector<size_t>& N) {

    // IC's
    T h0 = Height;
    T T0 = Angle;
    T tau = 0;
    int_ll x0 = 0;
    int bid = -1;
    bool sym = false;

    size_t j = 0;
    T TAU = 0;

    Scatterer<T> S(d, h0, T0, x0, bid);	// Evolves the map
    ScatterTraj<T> ST(N.size());

    for (size_t i = 0; i < N.back(); ++i) {
	if (N[j] == i) {
	    ST.Update(h0, T0, TAU, x0, bid); // 
	    j++;
	}
	S.Evolve_(h0, T0, x0, tau, bid, sym);
	TAU += tau;
	sym = false; // MUST RESET!!!!
    }
    ST.Update(h0, T0, TAU, x0, bid); // 
    return ST;
}

template <typename T>
int whichBranch(const T& Height, const T& Angle, const T& d) {
    for (size_t i = 0; i <= 19; ++i) {
	if (isPointInBranch(Height, Angle, {(int)i}, d))
	    return i;
    }
    std::cout << "@whichBranch(h,a,d) did not find a solution. This should be impossible?" << std::endl;
    return -1;
}

template <typename T>
ScatterVar<T> getPersistenceState(const T& Height, const T& Angle, const T& d, const std::vector<int>& IDs) {

    // IC's
    T h0 = Height;
    T T0 = Angle;
    T tau = 0;
    int_ll x0 = 0;
    int bid = -1;
    bool sym = false;

    size_t j = 0;
    T TAU = 0;

    Scatterer<T> S(d, h0, T0, x0, bid);	// Evolves the map
    ScatterVar<T> StatePrior(h0, T0, TAU, x0, bid);
    // Evolve whilst we are in a particular region or set of regions
    do {
	StatePrior = ScatterVar<T>(h0, T0, TAU, x0, bid);
	S.Evolve_(h0, T0, x0, tau, bid, sym);
	//std::cout << S << std::endl;
	//std::cout << bid << std::endl;
	TAU += tau;
	sym = false; // MUST RESET!!!!
    } while (In<int>(bid, IDs));

    //StatePrior.Print();
    return StatePrior;
}

template <typename T>
ScatterVar<T> JumpOnce(const T& Height, const T& Angle, const T& d) {

    // IC's
    T h0 = Height;
    T T0 = Angle;
    T tau = 0;
    int_ll x0 = 0;
    int bid = -1;
    bool sym = false;
    T TAU = 0;

    Scatterer<T> S(d, h0, T0, x0, bid);	// Evolves the map
    S.Evolve_(h0, T0, x0, tau, bid, sym);
    TAU += tau;

    return ScatterVar<T>(h0, T0, TAU, x0, bid);
}


#endif


