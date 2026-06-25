#ifndef SCATTER_H_INCLUDED
#define SCATTER_H_INCLUDED

const float_ldub PI = boost::math::constants::pi<float_ldub>();
const float_ldub P2 = PI/(float_ldub)2; 
const float_ldub A0 = 3*P2;
const float_ldub A20 = P2;

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

    

    std::vector<int> BranchID;      // A number between 0 - 19 indicating which 
    std::vector<int_ll> Position;	    // Lifted position of particle
    std::vector<T> Height, Angle, Tau;	    // Entrance heights, angles and jump lengths (cumulative)
    size_t idx;
};


// # Gamma 0 angles
//template <typename T> 
//T a0(const T& d, const T& h0)  { return (-PI/(T)2 + 2*PI); }
template <typename T> 
T a1(const T& d, const T& h0)  { return (-atan(h0/d) + 2*PI); }
template <typename T> 
T a2(const T& d, const T& h0)  { return (-atan((h0+d)/(T)2) + 2*PI); }
template <typename T> 
T a3(const T& h0)  { return (-atan(h0) + 2*PI); }
template <typename T> 
T a4(const T& h0)  { return (-atan(2*h0/(T)3) + 2*PI); }
template <typename T> 
T a5(const T& h0)  { return (-atan(h0/(T)2) + 2*PI); }   //# Note: Crossover between neg/pos angles
template <typename T> 
T a6(const T& d, const T& h0)  { return (+atan((d-h0)/(T)2)); }
template <typename T> 
T a7(const T& d, const T& h0)  { return (+atan(2*(d-h0)/(T)3)); }
template <typename T> 
T a8(const T& d, const T& h0)  { return (+atan(d-h0)); }
template <typename T> 
T a9(const T& d, const T& h0)  { return (+atan((2*d-h0)/(T)2)); }
template <typename T> 
T a10(const T& d, const T& h0) { return (+atan(2*(2*d-h0)/(T)3)); }
template <typename T> 
T a11(const T& d, const T& h0) { return (+atan(2*(3*d-h0)/(T)3)); } // # Note: angle doesn't coincide with the lower bound for a first bounce on Gamma2
										     //
// # Gamma 2 angles
template <typename T> 
T a12(const T& h0) { return atan(1 - 2*h0); }
template <typename T> 
T a13(const T& d, const T& h0) { return atan(1 - (2*h0/(1+2*d)) ); }
template <typename T> 
T a14(const T& d, const T& h0) { return atan(1 - h0/(d + 1)); }
template <typename T> 
T a15(const T& d, const T& h0) { return atan(1 + d - h0); }
template <typename T> 
T a16(const T& d, const T& h0) { return atan(1 + 2*(d - h0)); }

// # Gamma 3 angles
template <typename T> 
T a17(const T& d, const T& h0) { return atan(1 + 2*(d - h0)); }
template <typename T> 
T a18(const T& d, const T& h0) { return atan(2*d + 1 - h0); }
template <typename T> 
T a19(const T& d, const T& h0) { return atan(4*d + 1 - 2*h0); }
//template <typename T> 
//T a20(const T& d, const T& h0) { return PI/(T)2; }
template <typename T> 
T phi(const T& d, const T& h0) { return (-atan(2*(h0 + d)/(T)5) + 2*PI); }

template <typename T> 
int_ll kappa(const T& d, const T& h0, const T& t0) { 
    return ceil((d - h0 - tan(t0)) / (d*(tan(T(t0) - (T)1))));
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

	// @brief Iterates the scattering map one time
	void Evolve(T& h0, T& t0, int_ll& x0, T& tau, int& bid, bool& sym) {

	    if (h0 > d || h0 < 0) {
		std::cout << "h0 is out of bounds." << std::endl;
	    }

	    if ((t0 > 3*P2 && t0 < 2*PI) || (t0 >= 0 && t0 < atan(1- 2*h0))) { // BOTTOM LEFT
		Sg0(h0, t0, tau, bid);
		UpdatePosition(t0, x0, sym);
		nIterations++;
	    } else if (t0 > atan(1 - 2*h0) && t0 < atan(1 + 2*d - 2*h0)) {	// TOP RIGHT
		Sg2(h0, t0, tau, bid);
		UpdatePosition(t0, x0, sym);
		nIterations++;
	    } else if (t0 > atan(1 + 2*d - 2*h0) && t0 < P2) {		   // TOP LEFT
		Sg3(h0, t0, tau, bid);
		UpdatePosition(t0, x0, sym);
		nIterations++;
	    } else {
		sym = true;	
		if (t0 > PI && t0 < 3*P2) {
		    T temp = 3*PI - t0; // Reference variable
		    Evolve(h0, temp, x0, tau, bid, sym);
		    t0 = 3*PI - temp; // Revert 
		} else if (t0 > P2 && t0 < PI) {
		    T temp = PI - t0;
		    Evolve(h0, temp, x0, tau, bid, sym);
		    t0 = PI - temp;
		} else {
		    //t0 = fmod(t0, 2*PI); // Check this is the same as R's %% operation
		    //Evolve(h0, t0, x0, tau, sym);
		    std::cout << "-----------" << std::endl;
		    std::cout << nIterations << std::endl;
		    std::cout << "Height =   " << h0  << std::endl;
		    std::cout << "Angle =    " << t0 * 180/PI << std::endl;
		    std::cout << "Position = " << x0 << std::endl;
		    std::cout << "Branch =   " << bid << std::endl;
		    std::cout << "-----------" << std::endl;
		}

		if (t0 < 0)
		    t0 = fmod(t0 + 2*PI, 2*PI); // Check this is the same as R's %% operation
		else
		    t0 = fmod(t0, 2*PI);

		// When exiting the above Evolve(), the exit angle is returned to the
		// right half circle (positive horizontal velocity). Thus updating
		// vars here mean that we only record the incomplete scattering map.
		//Assign(h0, t0, x0, tau, bid);
	    }
	    // Update local vars for printing etc
	    Assign(h0, t0, x0, tau, bid);
	}

	void Update(T& h0, T& t0, int_ll& x0, T& tau, int& bid, bool& sym)  {

	    if ((t0 > 3*P2 && t0 < 2*PI) || (t0 >= 0 && t0 < atan(1- 2*h0))) { // BOTTOM LEFT
		Sg0(h0, t0, tau, bid);
		UpdatePosition(t0, x0, sym);
	    } else if (t0 > atan(1 - 2*h0) && t0 < atan(1 + 2*d - 2*h0)) {     // TOP RIGHT
		Sg2(h0, t0, tau, bid);
		UpdatePosition(t0, x0, sym);
	    } else {							       // TOP LEFT
		Sg3(h0, t0, tau, bid);
		UpdatePosition(t0, x0, sym);
	    } 
	    Assign(h0, t0, x0, tau, bid);
	}


	// @brief Iterates the scattering map one time WITHOUT RECURSION
	void Evolve_(T& h0, T& t0, int_ll& x0, T& tau, int& bid, bool& sym) {

	    if (h0 > d || h0 < 0) {
		std::cout << "h0 is out of bounds." << std::endl;
	    }

	    if ((t0 > 3*P2 && t0 < 2*PI) || (t0 >= 0 && t0 < P2)) { // Rightward moving
		Update(h0, t0, x0, tau, bid, sym);
		nIterations++;
	    } else {
		sym = true;	
		if (t0 > PI && t0 < (3*P2)) {
		    T temp = 3*PI - t0; // Reference variable
		    Update(h0, temp, x0, tau, bid, sym);
		    t0 = 3*PI - temp; // Revert 
		} else if (t0 > P2 && t0 <= PI) {
		    T temp = PI - t0;
		    Update(h0, temp, x0, tau, bid, sym);
		    t0 = PI - temp;
		} else {
		    std::cout << "-----------" << std::endl;
		    std::cout << nIterations << std::endl;
		    std::cout << "Height =   " << h0  << std::endl;
		    std::cout << "Angle =    " << t0 * 180/PI << std::endl;
		    std::cout << "Position = " << x0 << std::endl;
		    std::cout << "Branch =   " << bid << std::endl;
		    std::cout << "-----------" << std::endl;
		}

		if (t0 < 0) {
		    //std::cout << BranchID << std::endl;
		    t0 = fmod(t0 + 2*PI, 2*PI); // Check this is the same as R's %% operation
		} else {
		    t0 = fmod(t0, 2*PI);
		}
	    }
	}

	// @brief Returns index of the Sg0 branch
	void Sg0(T& h0, T& t0, T& tau, int& bid) {
	    if (t0 > A0 && t0 < a1(d,h0)) {
		tau = - h0 / sin(t0);
		h0 = -h0 * cot(t0);
		t0 = (5*P2) - t0;
		bid = 0;
	    } else if ((a2(d,h0) > a1(d,h0)) && (t0 > a1(d,h0) && t0 < a2(d,h0))) {
		tau = d*(alpha(d,h0,t0)+3) / cos(t0);
		h0 = h0 + d*(alpha(d,h0,t0) + 3)*tan(t0) + d*(alpha(d,h0,t0) + 1);
		t0 = t0;
		bid = 1;
	    } else if (((h0 > 1/(T)3) && (t0 > a2(d,h0) && t0 < a3(h0))) ||
		    (((h0 < 1/(T)3) && (a2(d,h0) >= a1(d,h0))) && (t0 > a2(d,h0) && t0 < phi(d,h0))) ||
		    (((h0 < 1/(T)3) && (a2(d,h0) < a1(d,h0))) && (t0 > a1(d,h0) && t0 < phi(d,h0)))) {
		tau = -(h0 + d) / sin(t0);
		h0 = (h0+d)*cot(t0) + 5/(T)2;
		t0 = (7*P2) - t0;
		bid = 2;
	    } else if ((((h0 < 1/(T)3) && (phi(d,h0) >= a1(d,h0))) && (t0 > phi(d,h0) && t0 < a3(h0))) ||
		    (((h0 < 1/(T)3) && (phi(d,h0) < a1(d,h0))) && (t0 > a1(d,h0) && t0 < a3(h0)))) {
		tau = 5*d/cos(t0);
		h0 = h0 + d + 5*tan(t0)/(T)2;
		t0 = t0;
		bid = 3;
	    } else if (t0 > a3(h0) && t0 < a4(h0)) {
		tau = 3*d/cos(t0);
		h0 = h0+(1+d)*tan(t0) + d;
		t0 = t0;
		bid = 4;
	    } else if (t0 > a4(h0) && t0 < a5(h0)) {
		tau = -h0/sin(t0);
		h0 = h0*cot(t0) + 4*d;
		t0 = (7*P2) - t0;
		bid = 5;
	    } else if ((t0 > a5(h0) && t0 < 2*PI) || (t0 >= 0 && t0 < a6(d,h0))) {
		tau = 4*d/cos(t0);
		h0 = h0 + 4*d*tan(t0);
		t0 = t0;
		bid = 6;
	    } else if (t0 > a6(d,h0)  && t0 < a7(d,h0)) {
		tau = 5*d/cos(t0);
		h0 = h0 + (2 + d)*tan(t0) - d;
		t0 = t0;
		bid = 7;
	    } else if (t0 > a7(d,h0)  && t0 < a8(d,h0)) {
		tau = 4*d/cos(t0);
		h0 = -h0 + 2*(d - tan(t0));
		t0 = t0 + PI;
		bid = 8;
	    } else if (t0 > a8(d,h0)  && t0 < a9(d,h0)) {
		tau = 4*d/cos(t0);
		h0 = -h0 + 1 - 2*tan(t0);
		t0 = t0 + PI;
		bid = 9;
	    } else if (t0 > a9(d,h0)  && t0 < a10(d,h0)) {
		tau = (2*d-h0)/sin(t0);
		h0 = 4*d - (1-h0)/tan(t0);
		t0 = (3*P2) - t0;
		bid = 10;
	    } else if (t0 > a10(d,h0) && t0 < a11(d,h0)) {
		tau = 3*d/cos(t0);
		h0 = -h0 + 3*d - (1 + d)*tan(t0);
		t0 = t0 + PI;
		bid = 11;
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
	    if (t0 > a12(h0) && t0 < a13(d,h0)) {
		    tau = (1+kappa(d,h0,t0)*d)/cos(t0);
		    h0 = h0 + d*kappa(d,h0,t0) * (tan(t0) - 1) + tan(t0);
		    t0 = t0;
		    bid = 12;
	    } else if (t0 > a13(d,h0) && t0 < a14(d,h0)) {
		    tau = 3*d/cos(t0);
		    h0 = -h0 + (1 - tan(t0))*(1 + d);
		    t0 = t0 + PI;
		    bid = 13;
	    } else if (t0 > a14(d,h0) && t0 < a15(d,h0)) {
		    tau = (3*d - h0)/sin(t0);
		    h0 = cot(t0) * (h0 + (d + 1)*(tan(t0) - 1));
		    t0 = -t0 + 3*P2;
		    bid = 14;
	    } else if (t0 > a15(d,h0) && t0 < a16(d,h0)) {
		    tau = 2*d/cos(t0);
		    h0 = -h0 + 2*d + 1 - tan(t0);
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
		    tau = 2*d/cos(t0);
		    h0 = 2*d + 1 - tan(t0) - h0;
		    t0 = t0 + PI;
		    bid = 16;
	    } else if (t0 > a18(d,h0) && t0 < a19(d,h0)) {
		    tau = (4*d-h0)/sin(t0);
		    h0 = cot(t0)*(h0 - 1 - 2*d + tan(t0));
		    t0 = (3*P2) - t0;
		    bid = 17;
	    } else if (t0 > a19(d,h0) && t0 < A20) { // Upper bound is a20
		    tau = (d*(gamma(d,h0,t0) + 3) - h0)/sin(t0);
		    h0 = cot(t0) * (-h0 + 1 + d*gamma(d,h0,t0)*(1-tan(t0)) + d*(1 + tan(t0)));
		    t0 = fmod((5*P2) - t0, 2*PI);
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

	// Given a (h0, t0) check that only one condition is satisfied each time
	void isUnique(const T& h0, const T& t0) {
	    int Sg0_counter = 0;
	    int Sg2_counter = 0;
	    int Sg3_counter = 0;

	    if ((t0 > 3*P2 && t0 < 2*PI) || (t0 >= 0 && t0 < atan(1- 2*h0))) { // BOTTOM LEFT
		if (t0 > A0 && t0 < a1(d,h0)) {
		    Sg0_counter++;
		} 
		if ((a2(d,h0) > a1(d,h0)) && (t0 > a1(d,h0) && t0 < a2(d,h0))) {
		    Sg0_counter++;
		} 
		if (((h0 > 1/(T)3) && (t0 > a2(d,h0) && t0 < a3(h0))) ||
			(((h0 < 1/(T)3) && (a2(d,h0) >= a1(d,h0))) && (t0 > a2(d,h0) && t0 < phi(d,h0))) ||
			(((h0 < 1/(T)3) && (a2(d,h0) < a1(d,h0))) && (t0 > a1(d,h0) && t0 < phi(d,h0)))) {
		    Sg0_counter++;
		} 
		if ((((h0 < 1/(T)3) && (phi(d,h0) >= a1(d,h0))) && (t0 > phi(d,h0) && t0 < a3(h0))) ||
			(((h0 < 1/(T)3) && (phi(d,h0) < a1(d,h0))) && (t0 > a1(d,h0) && t0 < a3(h0)))) {
		    Sg0_counter++;
		} 
		if (t0 > a3(h0) && t0 < a4(h0)) {
		    Sg0_counter++;
		} 
		if (t0 > a4(h0) && t0 < a5(h0)) {
		    Sg0_counter++;
		} 
		if ((t0 > a5(h0) && t0 < 2*PI) || (t0 >= 0 && t0 < a6(d,h0))) {
		    Sg0_counter++;
		} 
		if (t0 > a6(d,h0)  && t0 < a7(d,h0)) {
		    Sg0_counter++;
		} 
		if (t0 > a7(d,h0)  && t0 < a8(d,h0)) {
		    Sg0_counter++;
		} 
		if (t0 > a8(d,h0)  && t0 < a9(d,h0)) {
		    Sg0_counter++;
		} 
		if (t0 > a9(d,h0)  && t0 < a10(d,h0)) {
		    Sg0_counter++;
		} 
		if (t0 > a10(d,h0) && t0 < a11(d,h0)) {
		    Sg0_counter++;
		}

		if (Sg0_counter == 0)
		    Sg0_counter = -1;

	    } else if (t0 > atan(1 - 2*h0) && t0 < atan(1 + 2*d - 2*h0)) {     // TOP RIGHT
		if (t0 > a12(h0) && t0 < a13(d,h0)) {
		    Sg2_counter++;
		} 
		if (t0 > a13(d,h0) && t0 < a14(d,h0)) {
		    Sg2_counter++;
		} 
		if (t0 > a14(d,h0) && t0 < a15(d,h0)) {
		    Sg2_counter++;
		} 
		if (t0 > a15(d,h0) && t0 < a16(d,h0)) {
		    Sg2_counter++;
		}

		if (Sg2_counter == 0)
		    Sg2_counter = -1;
	    } else {							       // TOP LEFT
		if (t0 > a17(d,h0) && t0 < a18(d,h0)) {
		    Sg3_counter++;
		} 
		if (t0 > a18(d,h0) && t0 < a19(d,h0)) {
		    Sg3_counter++;
		} 
		if (t0 > a19(d,h0) && t0 < A20) { // Upper bound is a20
		    Sg3_counter++;
		} 
		
		if (Sg3_counter == 0)
		    Sg3_counter = -1;
	    }

	    // All 0
	    if (Sg0_counter == 0 && Sg2_counter == 0 && Sg3_counter == 0) {
		std::cout << "No interval was valid." << std::endl;
		std::cout << "h0 = " << h0 << ", theta0 = " << t0 *180/PI << std::endl;
	    }

	    // Any > 1
	    if (Sg0_counter > 1 || Sg2_counter > 1 || Sg3_counter > 1) {
		std::cout << "An interval was counted more than once." << std::endl;
		std::cout << "h0 = " << h0 << ", theta0 = " << t0 *180/PI << std::endl;
	    }

	    // Any > 1
	    if ((Sg0_counter + Sg2_counter + Sg3_counter) > 1) {
		std::cout << "Sum > 1" << std::endl;
		std::cout << "h0 = " << h0 << ", theta0 = " << t0 *180/PI << std::endl;
	    }

	    //std::cout << "Number of valid Sg0 branches = " << Sg0_counter << std::endl;
	    //std::cout << "Number of valid Sg2 branches = " << Sg2_counter << std::endl;
	    //std::cout << "Number of valid Sg3 branches = " << Sg3_counter << std::endl;
	}


	T getHeight() const { return H1; }
	T getAngle() const { return T1; }
	int_ll getPosition() const { return X1; }
	T getTime() const { return TAU; }
	size_t getBranchID() const { return BranchID; }
	size_t getIterations() const { return nIterations; }

    private:
	T d;
	size_t nIterations;
	size_t BranchID;    // A number between 0 - 19 indicating which 
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

#endif



// @Note. Initially, it was considered a good idea to check intervals independantly of updating the data. 
//	  Clearly, this means that every check gets done twice.
////// @brief For a first collision on Gamma0
//void Sg0(T& h0, T& t0, T& tau) {
//    BranchID = whichBranch_G0(d, h0, t0);
//    switch (BranchID) {
//	case 0: 
//	    tau = - h0 / sin(t0);
//	    h0 = -h0 * cot(t0);
//	    t0 = (5*P2) - t0;
//	    break;
//	case 1:
//	    tau = d*(alpha(d,h0,t0)+3) / cos(t0);
//	    h0 = h0 + d*(alpha(d,h0,t0) + 3)*tan(t0) + d*(alpha(d,h0,t0) + 1);
//	    t0 = t0;
//	    break;
//	case 2:
//	    tau = -(h0 + d) / sin(t0);
//	    h0 = (h0+d)*cot(t0) + 5/(T)2;
//	    t0 = (7*P2) - t0;
//	    break;
//	case 3:
//	    tau = 5*d/cos(t0);
//	    h0 = h0 + d + 5*tan(t0)/(T)2;
//	    t0 = t0;
//	    break;
//	case 4:
//	    tau = 3*d/cos(t0);
//	    h0 = h0+(1+d)*tan(t0) + d;
//	    t0 = t0;
//	    break;
//	case 5:
//	    tau = -h0/sin(t0);
//	    h0 = h0*cot(t0) + 4*d;
//	    t0 = (7*P2) - t0;
//	    break;
//	case 6:
//	    tau = 4*d/cos(t0);
//	    h0 = h0 + 4*d*tan(t0);
//	    t0 = t0;
//	    break;
//	case 7:
//	    tau = 5*d/cos(t0);
//	    h0 = h0 + (2 + d)*tan(t0) - d;
//	    t0 = t0;
//	    break;
//	case 8:
//	    tau = 4*d/cos(t0);
//	    h0 = -h0 + 2*(d - tan(t0));
//	    t0 = t0 + PI;
//	    break;
//	case 9:
//	    tau = 4*d/cos(t0);
//	    h0 = -h0 + 1 - 2*tan(t0);
//	    t0 = t0 + PI;
//	    break;
//	case 10:
//	    tau = (2*d-h0)/sin(t0);
//	    h0 = 4*d - (1-h0)/tan(t0);
//	    t0 = (3*P2) - t0;
//	    break;
//	case 11:
//	    tau = 3*d/cos(t0);
//	    h0 = -h0 + 3*d - (1 + d)*tan(t0);
//	    t0 = t0 + PI;
//	    break;
//	default: 
//	    std::cout << "Sg0 default branch. You will never see this message." << std::endl;
//    }
//}
//
//
//void Sg2(T& h0, T& t0, T& tau) {
//    BranchID = whichBranch_G2(d, h0, t0);
//    switch (BranchID) {
//	case 12: 
//	    tau = (1+kappa(d,h0,t0)*d)/cos(t0);
//	    h0 = h0 + d*kappa(d,h0,t0) * (tan(t0) - 1) + tan(t0);
//	    t0 = t0;
//	    break;
//	case 13:
//	    tau = 3*d/cos(t0);
//	    h0 = -h0 + (1 - tan(t0))*(1 + d);
//	    t0 = t0 + PI;
//	    break;
//	case 14:
//	    tau = (3*d - h0)/sin(t0);
//	    h0 = cot(t0) * (h0 + (d + 1)*(tan(t0) - 1));
//	    t0 = -t0 + 3*P2;
//	    break;
//	case 15:
//	    tau = 2*d/cos(t0);
//	    h0 = -h0 + 2*d + 1 - tan(t0);
//	    t0 = t0 + PI;
//	    break;
//	default: 
//	    std::cout << "Sg2 default branch. You will never see this message." << std::endl;
//    }
//}
//
//
//void Sg3(T& h0, T& t0, T& tau) {
//    BranchID = whichBranch_G3(d, h0, t0);
//    switch (BranchID) {
//	case 16:
//	    tau = 2*d/cos(t0);
//	    h0 = 2*d + 1 - tan(t0) - h0;
//	    t0 = t0 + PI;
//	    break;
//	case 17:
//	    tau = (4*d-h0)/sin(t0);
//	    h0 = cot(t0)*(h0 - 1 - 2*d + tan(t0));
//	    t0 = (3*P2) - t0;
//	    break;
//	case 18:
//	    tau = (d*(gamma(d,h0,t0) + 3) - h0)/sin(t0);
//	    h0 = cot(t0) * (-h0 + 1 + d*gamma(d,h0,t0)*(1-tan(t0)) + d*(1 + tan(t0)));
//	    t0 = fmod((5*P2) - t0, 2*PI);
//	    break;
//	default: 
//	    std::cout << "Sg3 default branch. You will never see this message." << std::endl;
//    }
//}
//
//// @brief Returns index of the Sg0 branch
//int whichBranch_G0(const T& d, const T& h0, const T& t0) {
//    if (t0 > A0 && t0 < a1(d,h0)) {
//	return 0;
//    } else if ((a2(d,h0) > a1(d,h0)) && (t0 > a1(d,h0) && t0 < a2(d,h0))) {
//	return 1;
//    } else if (((h0 > 1/(T)3) && (t0 > a2(d,h0) && t0 < a3(h0))) ||
//	    (((h0 < 1/(T)3) && (a2(d,h0) >= a1(d,h0))) && (t0 > a2(d,h0) && t0 < phi(d,h0))) ||
//	    (((h0 < 1/(T)3) && (a2(d,h0) < a1(d,h0))) && (t0 > a1(d,h0) && t0 < phi(d,h0)))) {
//	//			(( (h0 > 1/(T)3) && (t0 > a2(d,h0) && t0 < a3(h0))) 
//	//		    || ((h0 < 1/(T)3) && (t0 > a2(d,h0) && t0 < phi(d,h0)))) {
//	return 2;
//    } else if ((((h0 < 1/(T)3) && (phi(d,h0) >= a1(d,h0))) && (t0 > phi(d,h0) && t0 < a3(h0))) ||
//	    (((h0 < 1/(T)3) && (phi(d,h0) < a1(d,h0))) && (t0 > a1(d,h0) && t0 < a3(h0)))) {
//	// } else if ((h0 < 1/(T)3) && (t0 > phi(d,h0) && t0 < a3(h0))) {
//	return 3;
//    } else if (t0 > a3(h0) && t0 < a4(h0)) {
//	return 4;
//    } else if (t0 > a4(h0) && t0 < a5(h0)) {
//	return 5;
//    } else if ((t0 > a5(h0) && t0 < 2*PI) || (t0 >= 0 && t0 < a6(d,h0))) {
//	return 6;
//    } else if (t0 > a6(d,h0)  && t0 < a7(d,h0)) {
//	return 7;
//    } else if (t0 > a7(d,h0)  && t0 < a8(d,h0)) {
//	return 8;
//    } else if (t0 > a8(d,h0)  && t0 < a9(d,h0)) {
//	return 9;
//    } else if (t0 > a9(d,h0)  && t0 < a10(d,h0)) {
//	return 10;
//    } else if (t0 > a10(d,h0) && t0 < a11(d,h0)) {
//	return 11;
//    } else {
//	std::cout << "Impossble angle in Sg0" << std::endl;
//	std::cout << "h0 = " << h0 << std::endl;
//	std::cout << "t0 = " << t0*180/PI << std::endl;
//	return -1;
//    }
//    }
//
//    // @brief Returns index of the Sg2 branch
//    int whichBranch_G2(const T& d, const T& h0, const T& t0) {
//	if (t0 > a12(h0) && t0 < a13(d,h0)) {
//	    return 12;
//	} else if (t0 > a13(d,h0) && t0 < a14(d,h0)) {
//	    return 13;
//	} else if (t0 > a14(d,h0) && t0 < a15(d,h0)) {
//	    return 14;
//	} else if (t0 > a15(d,h0) && t0 < a16(d,h0)) {
//	    return 15;
//	} else {
//	    std::cout << "Impossible angle in Sg2" << std::endl;
//	    return -1;
//	}
//    }
//
//    // @brief Returns index of the Sg3 branch
//    int whichBranch_G3(const T& d, const T& h0, const T& t0) {
//	if (t0 > a17(d,h0) && t0 < a18(d,h0)) {
//	    return 16;
//	} else if (t0 > a18(d,h0) && t0 < a19(d,h0)) {
//	    return 17;
//	} else if (t0 > a19(d,h0) && t0 < A20) { // Upper bound is a20
//	    return 18;
//	} else {
//	    std::cout << "Impossible angle in Sg3" << std::endl;
//	    return -1;
//	}
//    }
//
