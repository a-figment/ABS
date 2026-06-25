#include <iostream>
#include <algorithm>
#include <chrono>
#include <thread>
#include "Billiard.h"
#include "Trajectory.h"
#include "OnlineStatistics.h"
#include "Misc.h"
#include "Helper.h"
#include "Ensemble.h"
#include "Simulation.h"
#include "DataKeys.h"
#include "Writer.h"
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


    Random<float_> Angle;
    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < Sim.nInitialConditions; ++i)
    {
	if (10 * i % Sim.nInitialConditions == 0 && i != 0)
	    std::cout << "Simulated " << i << " initial conditions." << std::endl;

	//Particle<T> P0 = Table.getInitialState(0, -1, Sim.TotalTime);
	//T th0 = Angle.getUniformRandom(0, 2*PI);
    	//Particle<T> P0(Point<float_>(0, 0.5), Point<float_>(cos(th0), sin(th0)), -1, -1, Sim.Time.back());

        //float_ th0 = Angle.getUniformRandom(-35*PI/180, 35*PI/180);
        float_ th0 = Angle.getUniformRandom(0, 2*PI);
        Particle<float_> P0(Point<float_>(-0.35, 0.2), Point<float_>(cos(th0), sin(th0)), -1, -1, 1000);

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


int main(int argc, char *argv[])
{
    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();

    // DEFINING DATA KEYS
    std::cout << "Arugments received: ";
    for (int i = 0; i < argc; ++i)
        std::cout << argv[i] << " ";
    std::cout << std::endl;
    std::cout << std::setprecision(std::numeric_limits<float_>::digits10) << std::showpoint;
    std::vector<int> KeyIndices;
    for (int i = 2; i < argc; ++i)
        KeyIndices.push_back(std::atoi(argv[i]));

    // DATA KEYS FOR PARALLEL CHANNELS
    DataKeys DK({"p", "a", "d", "r"}, {1000, 1000, 1000, 20 }, {"PARTITIONS", "ALPHA", "RATIO", "RADIUS" });
    std::string Key = DK.MakeKey(KeyIndices);
    DK.Print();

    // PARALLEL CHANNEL BILLIARD PARAMETERS
    std::vector<float_> alpha_rational = getRationalMultipleOfPi<float_>(7);
    std::vector<float_> alpha(alpha_rational.size() * 2);
    for (size_t i = 0; i < alpha.size(); ++i) {
        if (i < (alpha_rational.size()))
            alpha[i] = alpha_rational[i];
        else
            alpha[i] = (alpha_rational[i - alpha_rational.size()] * (PHI - (3 / ((float_)(5)))));
    }
    float_ dx = (1/((float_)2));
    float_ dy = dx * cos(alpha[KeyIndices[1]] / (float_)2) / sin(alpha[KeyIndices[1]] / (float_)2);
    //std::vector<float_> ratio = seq<float_>(1/(float_)4, 10, 40);
    
    // We may want to check ballistic fronts. So, make sure all the radii coincide with all the different ratios we want to consider
    std::vector<float_> ratio = seq<float_>(.4, 10, 40);
    std::vector<float_> porewidths(ratio.size(), 0); //= { - 1/((float_)2), (float_)1 };
    for (size_t i = 0; i < porewidths.size(); ++i) {
        porewidths[i] = ratio[i] * dy;
    }
    porewidths[0] = 0.4;
    // ====================================


    // Params
    float_ A = alpha[KeyIndices[1]];
    float_ beta = (PI - A)/2;
    float_ d = porewidths[KeyIndices[2]];


    // Distance from centre
    float_ w = 0.3;

    // Always determine radius from the central scatterer 
    float_ cy = dy - w;
    float_ cx = 0;
    //float_ radius = (1/((float_)2 * sin(A))) * (cx + dx + (cx + dx)*cos(A) - cy*sin(A));
    //float_ radius = (1/(float_)2) * (1 + 2*cx)*cos(A/(float_)2) - cy*sin(A/(float_)2); // Same as above but dx = 1/2 substituted
    float_ radius = (1/(float_)2) * cos(A/(float_)2) - cy*sin(A/(float_)2); // Same as above but dx = 1/2 substituted

    std::cout << "alpha = " << A * 180/PI << ", d = " << d << " radius = " << radius << std::endl;
    
    // Centers of each circle/arc seg
    std::vector<float_> xC = { -dx, 0, dx, dx, 0, -dx };
    std::vector<float_> yC = { w, dy-w, w, d+w, d+dy-w, d+w };

    // LEFT BOTTOM ARC SEG -- CHECKED
    Point<float_> LB1(xC[0] - radius*cos(beta), yC[0] - radius*sin(beta)); 
    ////Point<float_> LB1(-dx, yC[0] - radius); 
    Point<float_> LB2(xC[0] + radius*cos(beta), yC[0] - radius*sin(beta)); 
    //ArcSegment<float_> LeftBotArc(Point<float_>(xC[0], yC[0]), LB1, LB2);

    // MIDDLE BOTTOM SCATTER
    Circle<float_> MidBotScat(radius, Point<float_>(xC[1], yC[1]));
    Point<float_> MB1(xC[1] - radius*cos(beta), yC[1] + radius*sin(beta)); // cos(beta+pi) = -cos(beta), sin(beta+pi) = -sin(beta)
    Point<float_> MB2(xC[1] + radius*cos(beta), MB1.y); // cos(beta+pi) = -cos(beta), sin(beta+pi) = -sin(beta)


    // RIGHT BOTTOM ARC SEG
    Point<float_> RB1(LB1.x + 2*dx, LB1.y); 
    Point<float_> RB2(LB2.x + 2*dx, LB2.y);
    ArcSegment<float_> RightBotArc(Point<float_>(xC[2], yC[2]), RB1, RB2);
    
    //Point<float_> RB1(xC[0] - radius*cos(beta) + 2*dx, LB2.y); 
    //Point<float_> RB2(dx, yC[3] - radius);
    //ArcSegment<float_> RightBotArc(Point<float_>(xC[2], yC[2]), RB1, RB2);

    // Line segments of the channel are now floating. 
    LineSegment<float_> Gamma0(LB2, MB1);
    LineSegment<float_> Gamma1(MB2, RB1);

    // Periodic boundary No. 1
    //LineSegment<float_> Gamma2(Point<float_>(dx, yC[2]-radius), Point<float_>(dx, yC[3] - radius));

    // RIGHT TOP CIRCLE
    Circle<float_> RightTopScat(radius, Point<float_>(xC[3], yC[3]));
    Point<float_> RT1(xC[3] - radius*cos(beta), yC[3] - radius*sin(beta)); // Only left intersection needed
    Point<float_> RT2(xC[3] + radius*cos(beta), yC[3] - radius*sin(beta)); 

    // MIDDLE TOP ARC SEG
    Point<float_> MT1(MB2.x, MB2.y + d);
    Point<float_> MT2(MB1.x, MB1.y + d);
    ArcSegment<float_> MidTopArc(Point<float_>(xC[4], yC[4]), MT1, MT2);

    LineSegment<float_> Gamma3(RT1, MT1);

    // LEFT TOP CIRCLE
    Circle<float_> LeftTopScat(radius, Point<float_>(xC[5], yC[5]));
    Point<float_> LTS(LB2.x, LB2.y + d);
    Point<float_> LTS2(LB1.x, LB1.y + d);

    LineSegment<float_> Gamma4(MT2, LTS);

    // Periodic boundary No. 2
    //LineSegment<float_> Gamma5(Point<float_>(-dx, yC[5]-radius), Point<float_>(-dx, w-radius));

    // Shift the absorbing boundaries
    LineSegment<float_> Gamma2(RB2, RT2); 
    LineSegment<float_> Gamma5(Point<float_>(RT2.x - 2*dx, RT2.y), Point<float_>(RB2.x - 2*dx, RB2.y));
    std::cout << Gamma2 << " " << Gamma5 << std::endl;

    Polygon<float_> Skeleton( {Gamma0, Gamma1, Gamma2, Gamma3, Gamma4, Gamma5}, {2, 5}); 


    //std::cout << LB2 << " " << MB1 << " " << MB2 << " " << RB1 << " " << RB2 << std::endl;
    
    // Writing all the data
    std::string DataDir = "Data/" + Key;
    Skeleton.Write(DataDir + "BoundaryPoints.dat");
    Logger<float_> Log(DataDir + "BoundaryPoints.dat", true);
    for (size_t i = 0; i < xC.size(); ++i) 
	Log.Stream << xC[i] << " " << yC[i] << std::endl;
    Log.Stream << radius << " " << radius << std::endl;
    std::cout << Skeleton << std::endl;

    Logger<float_> Log_Circs(DataDir + "3PtCircles.dat", true);

    Log_Circs.Stream << xC[2] - 1 << " " << yC[2] << std::endl;
    Log_Circs.Stream << xC[0] - 1 << " " << yC[0] << std::endl;
    Log_Circs.Stream << xC[2] - 1 << " " << yC[2] << std::endl;

    // Lower segments
    std::vector<std::unique_ptr<Obstacle<float_>>> Obstacles(5);
    Obstacles[0] = std::make_unique<Polygon<float_>>(Skeleton);
    //Obstacles[1] = std::make_unique<ArcSegment<float_>>(LeftBotArc);
    Obstacles[1] = std::make_unique<Circle<float_>>(MidBotScat);
    Obstacles[2] = std::make_unique<ArcSegment<float_>>(RightBotArc);
    Obstacles[3] = std::make_unique<Circle<float_>>(RightTopScat);
    Obstacles[4] = std::make_unique<ArcSegment<float_>>(MidTopArc);
    //Obstacles[6] = std::make_unique<Circle<float_>>(LeftTopScat);
    Billiard<float_> Table(Obstacles);

    //TODO: Other angles not working. Plotting fine, but boundary points are incorrect (tangent intersections).


    //Random<float_> Angle;
    //float_ th0 = Angle.getUniformRandom(0, 2*PI);
    //Particle<float_> P(Point<float_>(-0.05, 0.6), Point<float_>(cos(th0), sin(th0)), -1, -1, 1000);

    //Trajectory<float_> dTraj = Table.getTrajectory(P, 20);
    //Trajectory<float_> dTraj = Table.getTrajectory(P, Time);
    //std::vector<float_> x = dTraj.getLiftedPositions_x();
    //std::vector<float_> y = dTraj.getLiftedPositions_y();
    //std::vector<float_> fps = dTraj.getLiftedFreePaths();
    //std::vector<float_> t(fps.size(), 0);
    //std::partial_sum(fps.begin(), fps.end(), t.begin());

    //std::string DataDir = "Data/" + Key;
    //WriteColVector<float_>(DataDir + "t.dat", t);
    //WriteColVector<float_>(DataDir + "x.dat", x);
    //WriteColVector<float_>(DataDir + "y.dat", y);


    // TODO: Need to generate uniform positions on the cell

    // Animations
    //Random<float_> Angle;
    ////std::vector<float_> th0 = seq<float_>(-31*PI/180, -25*PI/180, 2000);
    ////std::vector<float_> th0 = { - 29.27713 * PI/180 };
    std::vector<float_> Time = seq<float_>(0, 2, 201);
    //for (size_t i = 0; i < 10000; ++i) {
    //    //float_ th0 = Angle.getUniformRandom(-29*PI/180, -28*PI/180);
    //    //Particle<float_> P(Point<float_>(0, 0.6), Point<float_>(cos(th0[i]), sin(th0[i])), -1, -1, 1000);
    //    
    //    float_ th0 = Angle.getUniformRandom(-35*PI/180, 35*PI/180);
    //    Particle<float_> P(Point<float_>(-0.35, 0.2), Point<float_>(cos(th0), sin(th0)), -1, -1, 1000);
    //    //Particle<float_> P(Point<float_>(MT2.x-0.1, MT2.y-0.15), Point<float_>(cos(PI/4), sin(PI/4)), -1, -1, 1000);

    //    //std::cout << th0 << std::endl;
    //    //std::cout << P.Pos << std::endl;
    //    Trajectory<float_> dTraj = Table.getTrajectory(P, Time);
    //    //Trajectory<float_> dTraj = Table.getTrajectory(P, 100);

    //    std::vector<float_> x = dTraj.getLiftedPositions_x();
    //    std::vector<float_> y = dTraj.getLiftedPositions_y();
    //    //std::vector<float_> x = dTraj.getPositions_x();
    //    //std::vector<float_> y = dTraj.getPositions_y();

    //    WriteVector<float_>(DataDir + "x.dat", x);
    //    WriteVector<float_>(DataDir + "y.dat", y);
    //    //std::vector<float_> fps = dTraj.getLiftedFreePaths();
    //    //std::vector<float_> t(fps.size(), 0);
    //    //std::partial_sum(fps.begin(), fps.end(), t.begin());
    //    //WriteVector<float_>(DataDir + "t.dat", t);
    //}

    //Writer<float_> w_X(DataDir + "x.dat", true);
    //Writer<float_> w_Y(DataDir + "y.dat", true);
    //for (size_t i = 0; i < 5000; ++i) {

    //    Particle<float_> P = Table.getIC_Rectangle(0, -1, Time.back(), Point<float_>(-0.4, 0.25), Point<float_>(-0.3, 0.35));
    //    P.Vel = Point<float_>(cos(20*PI/180), sin(20*PI/180));
    //	Trajectory<float_> tTraj = Table.getTrajectory(P, Time);

    //    w_X.WriteRowVector(tTraj.getLiftedPositions_x());
    //    w_Y.WriteRowVector(tTraj.getLiftedPositions_y());
    //}

    Writer<float_> w_X(DataDir + "x1.dat", true);
    Writer<float_> w_Y(DataDir + "y1.dat", true);
    Particle<float_> P;
    //P.Vel = Point<float_>(cos(90*PI/180), sin(90*PI/180));
    ////P.Vel = Point<float_>(1, 0);
    //P.Pos = Point<float_>(.5, 0.33);

    float_ x0 = -0.2;
    float_ y0 = yC[1]+0.4;
    P.Pos = Point<float_>(x0, y0);
    float_ angleToVertex = atan((radius + yC[1] - y0)/(-x0));
    float_ epsil = -2*PI/(float_)180;
    float_ An = atan((radius + yC[1] - y0)/(-x0));
    P.Vel = Point<float_>(cos(An), sin(An));
    P.Vel = Point<float_>(cos(An+epsil), sin(An+epsil));

    Trajectory<float_> tTraj = Table.getTrajectory(P, 20);
    w_X.WriteRowVector(tTraj.getLiftedPositions_x());
    w_Y.WriteRowVector(tTraj.getLiftedPositions_y());

    //Writer<float_> w_X(DataDir + "x2.dat", true);
    //Writer<float_> w_Y(DataDir + "y2.dat", true);
    //Particle<float_> P;
    ////P.Vel = Point<float_>(cos(1*PI/180), sin(1*PI/180));
    //P.Vel = Point<float_>(cos(87*PI/180), sin(87*PI/180));
    ////P.Vel = Point<float_>(1, 0);
    //P.Pos = Point<float_>(.5, 0.33);
    //Trajectory<float_> tTraj = Table.getTrajectory(P, 30);
    //w_X.WriteRowVector(tTraj.getLiftedPositions_x());
    //w_Y.WriteRowVector(tTraj.getLiftedPositions_y());


    // PARAMETERS FOR SIMULATION OBJECT
    //std::vector<float_> T = { 1e4, 1e5, 1e6, 1e7 }; // Total time
    //std::vector<uint> nSteps = { 101, 1001 };  // Number of time steps
    //std::vector<float_> TIME = seq<float_>(0, T[1], nSteps[1]); // Time
    //std::vector<float_> ttg = { TIME[1], TIME[10], TIME[100], TIME.back() };	// Times to inspect an observable
    //int M = 100000;	      // Number of initial conditions
    //int nm = 2;		      // Number of moments
    //float_ mstep = (float_)1; // Moment step
    //
    ////std::string DataDir = "Data/" + Key;
    ////std::string DataDir = "/fred/oz055/jorchard_data/ABS/S-RoundedCorners/" + Key;
    //
    //Simulation<float_> Sim {
    //    .Time = TIME,
    //    .L = {0},
    //    .TimesToGet = ttg,
    //    .Symbols = {""},
    //    .DataDestination = DataDir,
    //    .RunTime = ".",
    //    .MomentStep = mstep,
    //    .TotalTime = TIME.back(),
    //    .nInitialConditions = M,
    //    .nCollisions = 0, // always reach max time before max coll
    //    .nMoments = nm};
    //Sim.Print();

    //getMoments_t<float_>(Table, Sim);

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    const std::string RunTime = "\033[1;32mSimulation run time: " + std::to_string(std::chrono::duration_cast<std::chrono::hours>(t1 - t0).count()) + " hours = " + std::to_string(std::chrono::duration_cast<std::chrono::minutes>(t1 - t0).count()) + " minutes = " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) + " milliseconds. \n\033[0m" + "\n";
    return 0;
}


    

