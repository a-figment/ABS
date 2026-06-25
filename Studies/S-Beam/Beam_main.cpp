#include <iostream>
#include <algorithm>
#include <chrono>
#include <thread>
#include "Billiard.h"
#include "Trajectory.h"
#include "OnlineStatistics.h"
#include "Misc.h"
#include "Ensemble.h"
#include "Simulation.h"
#include "DataKeys.h"
#include "Logger.h"

// < single character prefix, max number of instantiations >
typedef float_ldub float_;
const int nDigits = std::numeric_limits<float_>::digits10;
const float_ PHI = boost::math::constants::phi<float_>();
const float_ PI = boost::math::constants::pi<float_>();


// @brief Writes moments of displacement and displacements at selected times
template <typename T>
inline void getMoments_t(Billiard<T> &Table, Simulation<T> &Sim) 
{
    // Involves all initial conditions
    std::vector<Ensemble<T>> Moments(Sim.nMoments, Ensemble<T>(Sim.Time.size()));

    Writer<T> w_DX(Sim.DataDestination + "DX.dat", true);
    std::vector<size_t> Indices = getIndicesOfVectorValue(Sim.TimesToGet, Sim.Time);
    std::cout << "Collecting displacements at times t = " << index(Sim.Time, Indices) << " with indices = " << Indices << std::endl;

    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < Sim.nInitialConditions; ++i)
    {
	if (10 * i % Sim.nInitialConditions == 0 && i != 0)
	    std::cout << "Simulated " << i << " initial conditions." << std::endl;

	Particle<T> P0 = Table.getInitialState(0, -1, Sim.TotalTime);

	Trajectory<T> Traj = Table.getTrajectory(P0, Sim.Time);
        std::vector<T> xDis = Traj.getDisplacements_x();
	w_DX.WriteRowVector(xDis, Indices);

	std::vector<std::vector<T>> DXp = Traj.getDisplacements_x(Sim.MomentStep, Sim.nMoments);
	for (int j = 0; j < Sim.nMoments; ++j)
	    Moments[j].Push(DXp[j]);

    }
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";

    for (size_t i = 0; i < Moments.size(); ++i)
	Moments[i].Write_Append(Sim.DataDestination + "Moments.dat");
}


template <typename T>
void WriteInitialBeamDistribution(Billiard<T>& Table) {
    size_t nBeams = 100000;
    size_t nParticlesOnBeam = 50; 
    size_t factor = 150;
    for (size_t i = 0; i < nBeams; ++i) {
      Beam<T> b = Table.getInitialBeam();
      nParticlesOnBeam = std::round(b.Length * factor); // n propto beam length
      b.Write("Data/Beams.dat");
    }
}

template <typename T>
void WriteFirstJump(Billiard<T>& Table) {
    Writer<T> w_FJ("Data/FJ.dat", true);
    size_t nBeams = 1000000;
    for (size_t i = 0; i < nBeams; ++i) {
      //Beam<T> b = Table.getInitialBeam();
      Beam<T> b = Table.getInitialBeamWithPB();

      T xi = b.B.x - b.A.x;
      T tau = b.Length;

      b.Write("Data/Beams.dat");
      w_FJ.WriteRowVector({ tau, xi, b.RayAngle });
    }
}

template <typename T>
void EvolveOneBeam(Billiard<T>& Table, Simulation<T>& Sim) {
    size_t factor = 1000;

    // Get one beam
    Beam<T> B = Table.getInitialBeam();
    size_t nPointsOnBeam = std::round(B.Length * factor); // n propto beam length

    // Evolve it and write positions
    Point<T> Vel(cos(B.VelAngle), sin(B.VelAngle));
    for (size_t i = 1; i < nPointsOnBeam; ++i)
    {
        T t = (T)i / (T)nPointsOnBeam;
        Point<T> Pos = B.getPointOnBeam(t);
        Particle<T> P(Pos, Vel, 0, -1, Sim.TotalTime);

	// Discrete
        //Trajectory<T> Traj = Table.getTrajectory(P, Sim.nCollisions);
	//Traj.WritePositions_App(Sim.DataDestination + std::to_string(i) + "-dTraj.dat");

	// Continuous
        Trajectory<T> Traj = Table.getTrajectory(P, Sim.Time);
	Traj.WritePositions_App(Sim.DataDestination + std::to_string(i) + "-cTraj.dat");
    }
}


template <typename T>
void JumpCell(Billiard<T>& Table, Particle<T>& P) {
    int currCID = P.CellID.x;
    while (P.CellID.x == currCID) 
	Table.Update(P);
}

// @brief Fix an initial velocity, and evolve all particles with uniform positions
//	  until they jump from one elementary cell to another. Then, 
// @goals Visualise slicing on an array of channels
template <typename T>
void EvolveCloud(Billiard<T>& Table, Simulation<T>& Sim) {
    // Writers
    Writer<T> w_X(Sim.DataDestination + "x.dat", true);
    Writer<T> w_Y(Sim.DataDestination + "y.dat", true);

    // Geometry
    Point<T> diam = Table.getDiameter();

    // Fixed uniform velocity
    Random<T> R;
    T theta = R.getUniformRandom(static_cast<T>(0), static_cast<T>(2) * boost::math::constants::pi<T>());
    Point<T> Vel(cos(theta), sin(theta)); // Uniform random angle
    std::cout << theta << std::endl;

    // Generate uniform initial positions and fixed velocity particles
    std::vector<T> xPos(Sim.nInitialConditions, 0);
    std::vector<T> yPos(Sim.nInitialConditions, 0);
    std::vector<Particle<T>> Ps(Sim.nInitialConditions, Particle<T>());
    for (size_t i = 0; i < (size_t)Sim.nInitialConditions; ++i)
    {
	Particle<T> P = Table.getInitialState(0, -1, Sim.TotalTime);
	P.Vel = Vel;
	Ps[i] = P;
	xPos[i] = P.Pos.x;
	yPos[i] = P.Pos.y;
    }

    size_t nCellUpdates = 10; // Number of times all particles change cells
    for (size_t k = 1; k < nCellUpdates; k++) {
	if (10 * k % nCellUpdates == 0 && k != 0)
	    std::cout << "Simulated " << k << " cell hops." << std::endl;

	// Write positions after each updte
	w_X.WriteRowVector(xPos);
	w_Y.WriteRowVector(yPos);

	// Evolve each particle until it changes elementary cell
	for (size_t i = 0; i < Ps.size(); ++i) {
	    JumpCell(Table, Ps[i]);
	    xPos[i] = Ps[i].LPos(diam).x;
	    yPos[i] = Ps[i].Pos.y;
	    //Ps[i] = P;
	}
    }
}

std::vector<int> which_equal(const std::vector<int>& ints, const int& test_val) {
    std::vector<int> result;
    result.reserve(ints.size());
    for (size_t i = 0; i < ints.size(); ++i)
	if (ints[i] == test_val)
	    result.push_back(ints[i]);
    return result;
}

// @brief Fix an initial velocity, and evolve all particles with uniform positions
//	  until they jump from one elementary cell to another. Then, 
// @goals Visualise slicing on an array of channels
template <typename T>
void DiscreteMoments(Billiard<T>& Table, Simulation<T>& Sim, const size_t& TimeSteps) {

    // Geometry
    Point<T> diam = Table.getDiameter();

    // Moments of displacement
    std::vector<int> p = seq<int>(2, 2*Sim.nMoments, Sim.nMoments); // Even integer moments
    //std::vector<std::vector<T>> Moments(p.size(), std::vector<T>(TimeSteps,0));
    std::vector<std::vector<size_t>> Moments(p.size(), std::vector<size_t>(TimeSteps,0));

    // Fixed uniform velocity
    //Random<T> R;
    //T theta = R.getUniformRandom(static_cast<T>(0), static_cast<T>(2) * boost::math::constants::pi<T>());
    //Point<T> Vel(cos(theta), sin(theta)); // Uniform random angle
    //std::cout << theta << std::endl;

    // Generate uniform initial positions and fixed velocity particles
    std::vector<Particle<T>> Ps(Sim.nInitialConditions, Particle<T>());
    for (size_t i = 0; i < (size_t)Sim.nInitialConditions; ++i)
    {
	Particle<T> P = Table.getInitialState(0, -1, Sim.TotalTime);
	//P.Vel = Vel;
	Ps[i] = P;
    }

    Logger<T> w_HC(Sim.DataDestination + "HeadingCounts.dat", true);
    //std::vector<T> MSD(TimeSteps, 0);
							 
    for (size_t k = 1; k < TimeSteps; k++) {
	if (10 * k % TimeSteps == 0 && k != 0) 
	    std::cout << "Simulated " << k << " cell hops." << std::endl;

	// Track the cell id of each particle
	std::vector<int> rounded_xpos(Ps.size(), 0);
	std::vector<T> transfer_time(Ps.size(), 0);

	// Evolve each particle until it changes elementary cell
	for (size_t i = 0; i < Ps.size(); ++i) {
	    JumpCell(Table, Ps[i]);
	    rounded_xpos[i] = std::round(Ps[i].LPos(diam).x);
	    transfer_time[i] = Ps[i].TotalTime;
	}
	
	if (k == TimeSteps - 1) {
	    std::cout << transfer_time << std::endl;
	    std::cout << rounded_xpos << std::endl;
	}

	// Discrete MSD
	std::vector<size_t> Counts(2*TimeSteps-1, 0); // Number of particles in 
	for (int j = -k; j <= (int)k; j+=2) { // Time is distance
	    std::vector<int> we = which_equal(rounded_xpos, j); // Which particles in jth cell
	    size_t Aj = we.size();			        // Number of particles in jth cell 

	    for (size_t i = 0; i < p.size(); ++i)  
		Moments[i][k] += (((T)Aj / (T)Sim.nInitialConditions) * (T)pow(j, p[i])); 
	    Counts[j+TimeSteps-1] = Aj;
	}
	//w_HC.Stream << Counts << std::endl;
    }

    // Write moments of displacmeent
    Writer<size_t> w_Moments(Sim.DataDestination + "Moments.dat", false);
    for (size_t i = 0; i < Moments.size(); ++i) 
	w_Moments.WriteRowVector(Moments[i]);

}

int main(int argc, char *argv[])
{
    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    // ======================================= BOILER PLATE START =======================================
    std::cout << "Arugments received: ";
    for (int i = 0; i < argc; ++i)
        std::cout << argv[i] << " ";
    std::cout << std::endl;

    std::cout << std::setprecision(std::numeric_limits<float_>::digits10) << std::showpoint;
    std::vector<int> KeyIndices;
    for (int i = 2; i < argc; ++i)
        KeyIndices.push_back(std::atoi(argv[i]));
    DataKeys DK({"p", "b", "t", "d"}, {1000, 40, 40, 40}, {"PARTITIONS", "ALPHA 1 (BOT)", "ALPHA 2 (TOP)", "PORE WIDTH"});
    std::string Key = DK.MakeKey(KeyIndices);
    DK.Print();

    std::vector<std::unique_ptr<Obstacle<float_>>> Obstacles(1);
    // ======================= BOILER PLATE END ==========================
    
    // Some rational and irrational angles
    std::vector<float_> alpha_rational = getRationalMultipleOfPi<float_>(1);
    std::vector<float_> alpha(alpha_rational.size() * 2);
    for (size_t i = 0; i < alpha.size(); ++i) {
        if (i < (alpha_rational.size()))
            alpha[i] = alpha_rational[i];
        else
            alpha[i] = (alpha_rational[i - alpha_rational.size()] * (PHI - (3 / ((float_)(5)))));
    }
    std::vector<float_> ratio = seq<float_>(1/(float_)4, 10, 976); // Steps of 0.01. 0.25, 0.26,...

    // ======================= BILLIARD DEFINITION START ==========================
    // PARALLEL CHANNELS
    //std::cout << "Initialising parallel channel with Orchard parameters (superdiffusion)" << std::endl;
    //float_ dx = (1/((float_)2));
    //float_ dy = dx * cos(alpha[KeyIndices[1]] / (float_)2) / sin(alpha[KeyIndices[1]] / (float_)2);
    //std::vector<float_> porewidths(ratio.size(), 0); //= { - 1/((float_)2), (float_)1 };
    //for (size_t i = 0; i < porewidths.size(); ++i) 
    //    porewidths[i] = ratio[i] * dy;
    //Polygon<float_> Channel = getParallelChannel2<float_>(dy, (1 / (float_)2), porewidths[KeyIndices[2]]);


    // NON-PARALLEL CHANNELS (SANDERS PARAMS)
    //std::cout << "Initialising non-parallel channel with Sanders' parameters" << std::endl;
    //std::vector<float_> H = { 207 / ((float_)100), 127 / ((float_)100), 117 / ((float_)100), 107 / ((float_)100) };
    //float_ dyb = 45/((float_)100);
    //float_ dyt = 77/((float_)100);
    //Polygon<float_> Channel = getNonParallelChannel<float_>((float_)1, dyb, dyt, H[KeyIndices[1]] - dyt);

    // NON-PARALLEL CHANNELS (SUBDIFFUSION, RONDONI PARAMETERS)
    //std::cout << "Initialising non-parallel channel with Rondoni parameters (sub-diffusion)" << std::endl;
    //std::vector<float_> H = { (float_)2, 3/(float_)2, (float_)1, 3/(float_)4, 1/(float_)2 };
    //float_ dyt = dx;
    //std::vector<float_> dyb_on_dx = {(float_)1 + 1/(float_)100, 
    //    			     (float_)1 + 1/(float_)1000, 
    //    			     (float_)1 + 1/(float_)10000, 
    //    			     (float_)1 + 1/(float_)100000, 
    //    			     (float_)1 + 1/(float_)1000000, 
    //    			     (float_)1}; // = dyb / dx 
    //std::vector<float_> dyb(dyb_on_dx.size());
    //for (size_t i = 0; i < dyb.size(); ++i) 
    //    dyb[i] = dyb_on_dx[i] * dx;
    //Polygon<float_> Channel = getNonParallelChannel<float_>(1/(float_)2, dyb[KeyIndices[2]], dyt, H[KeyIndices[1]] - dyt);


    // CHANNEL ANTI-SYMMETRIC
    //T centreline = 1/(float_)2;
    //std::vector<float_> w = { 1/(float_)10, 3/(float_)10 };

    // CHANNEL (TWO-ANGLE DEFINITION)
    if ((KeyIndices[1] == 0 && KeyIndices[2] == 2) || (KeyIndices[1] == 2 && KeyIndices[2] == 0)) // Invalid pairs
	return 0;
    if ((KeyIndices[1] == 3 && KeyIndices[2] == 5) || (KeyIndices[1] == 5 && KeyIndices[2] == 3)) // Invalid pairs
	return 0;
    bool isParallel = (KeyIndices[1] == KeyIndices[2]);
    std::cout << KeyIndices[1] << " " << KeyIndices[2] << std::endl;
    Polygon<float_> Channel = getChannel<float_>(1/(float_)2, alpha[KeyIndices[1]], alpha[KeyIndices[2]], ratio[KeyIndices[3]], isParallel);


    Obstacles[0] = std::make_unique<Polygon<float_>>(Channel);
    Billiard<float_> Table(Obstacles);
    // ======================= BILLIARD DEFINITION END ==========================
    
    std::string DataDir = "Data/" + Key;
    //std::string DataDir = "/fred/oz055/jorchard_data/ABS/S-Beam/" + Key;

    std::vector<float_> TIME = seq<float_>(0, 1e8, 1001);
    //std::vector<float_> TIME = seq<float_>(1/((float_)100), 100,  5000);
    Simulation<float_> Sim;
    Sim = Simulation<float_> {
        .Time = TIME,
        .L = { 0 },
        .TimesToGet = { 0 },
        .Symbols = {""},
        .DataDestination = DataDir,
        .RunTime = ".",
        .MomentStep = 0,
        .TotalTime = TIME.back(),
        .nInitialConditions = 10,
        .nCollisions = 0,
        .nMoments = 10};

    //WriteInitialBeamDistribution<float_>(Table);
    //WriteFirstJump<float_>(Table);
    //EvolveOneBeam<float_>(Table, Sim);
    //EvolveCloud<float_>(Table, Sim);
    DiscreteMoments<float_>(Table, Sim, 100001);

    //float_ factor = 10000;
    //Beam<float_> b = Table.getInitialBeam();
    //int nP = std::round(b.Length * factor);
    //std::cout << "Beam 0: " << b << std::endl;
    //std::cout << "Beam with " << nP << " particles." << std::endl;
    //EvolveBeam<float_>(b, Table, Sim, nP, 0);


    // CONTINUOUS TIME MOMENTS FOR COMPARISON
    std::vector<float_> ttg = { TIME[10], TIME[100], TIME[1000] };	// Times to inspect an observable
    Simulation<float_> MomentSim {
        .Time = TIME,
        .L = {0},
        .TimesToGet = ttg,
        .Symbols = {""},
        .DataDestination = DataDir,
        .RunTime = ".",
        .MomentStep = (float_)1,
        .TotalTime = TIME.back(),
        .nInitialConditions = 100000,
        .nCollisions = 0, // always reach max time before max coll
        .nMoments = 10};
    MomentSim.Print();
    //getMoments_t<float_>(Table, MomentSim);

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";
    std::cout << RunTime << std::endl;
    return 0;
}

