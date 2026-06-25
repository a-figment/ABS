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

// < single character prefix, max number of instantiations >
typedef float_ldub float_;
const int nDigits = std::numeric_limits<float_>::digits10;

int main(int argc, char *argv[])
{
    std::cout << "Arugments received: ";
    for (int i = 0; i < argc; ++i) 
        std::cout << argv[i] << " ";
    std::cout << std::endl;

    // ======================= BOILER PLATE START ==========================
    std::cout << std::setprecision(std::numeric_limits<float_>::digits10) << std::showpoint;
    std::vector<int> KeyIndices;
    for (int i = 2; i < argc; ++i) 
        KeyIndices.push_back(std::atoi(argv[i]));
    DataKeys DK({"p", "b"}, {1000, 5}, {"Partitions", "Boundary parameter"});
    std::string Key = DK.MakeKey(KeyIndices);
    DK.Print();
    // ======================= BOILER PLATE END ==========================
    
    // ======================= BILLIARD DEFINITION START ==========================
    
    std::string DataDestination = "Data/" + Key;
    // Internal angle definitions
    //float_ PHI = boost::math::constants::phi<float_>();
    float_ PI = boost::math::constants::pi<float_>();

    

    // Barrier billard (test)
    //LineSegment<float_> LS1(Point<float_>(0,0), Point<float_>(1, 0));
    //LineSegment<float_> LS2(Point<float_>(1,0), Point<float_>(1, 1));
    //LineSegment<float_> LS3(Point<float_>(1,1), Point<float_>(0, 1));
    //LineSegment<float_> LS4(Point<float_>(0,1), Point<float_>(0, 0));
    //LineSegment<float_> LS5(Point<float_>(0.25,0.5), Point<float_>(0.75, 0.5));
    //std::vector<LineSegment<float_>> Lines = { LS1, LS2, LS3, LS4, LS5 };
    // ======================

    std::vector<LineSegment<float_>> Lines;// = { LS1, LS2, LS3, LS4, LS5 };
    Lines.reserve(1000);

    // Entrance and exit segments
    float_ side_len = 3; // Assume our maze is always square and defined on integer lattice
    float_ gap = 1; // Distance between wall
    LineSegment<float_> Start(Point<float_>(0, side_len - gap), Point<float_>(0, side_len));
    LineSegment<float_> End(Point<float_>(side_len, 0), Point<float_>(side_len, gap));


    // Read maze
    std::fstream MazeLines("lines.dat", std::ios::in);  
    if (MazeLines.is_open()) {   //checking whether the file is open
	std::string line;
	while(getline(MazeLines, line)) { //read data from file object and put it into string.

	    // String of ints to vector of ints
	    std::stringstream iss(line);
	    int number;
	    std::vector<int> Endpoints;
	    while (iss >> number) {
		Endpoints.push_back(number);
	    }
	    std::cout << Endpoints << std::endl;

	    // Create each line segment of the maze
	    LineSegment<float_> LS(Point<float_>(Endpoints[0], Endpoints[1]), 
				   Point<float_>(Endpoints[2], Endpoints[3]));
	    Lines.push_back(LS);
	}

	// Close off the start and end of the maze
	//LineSegment<float_> Start(Point<float_>(10, 100), Point<float_>(10, 110));
	//LineSegment<float_> End(Point<float_>(110, 10), Point<float_>(110, 20));

	//LineSegment<float_> Start(Point<float_>(0, 0), Point<float_>(0, 1000));
	//LineSegment<float_> End(Point<float_>(120, 0), Point<float_>(120, 1000));
	Lines.push_back(Start);
	Lines.push_back(End);
	MazeLines.close(); //close the file object.
    }


    FloatingPolygon<float_> Maze(Lines, { -1, -1 });
    std::cout << Maze << std::endl;

    std::vector<std::unique_ptr<Obstacle<float_>>> Obstacles(1);
    Obstacles[0] = std::make_unique<FloatingPolygon<float_>>(Maze);
    Billiard<float_> Table(Obstacles);

    //Particle<float_> P(Point<float_>(0,0), Point<float_>(cos(PI/4), sin(PI/4)), -1, 10);
    // How to make a floating line segment which you can collide with from both directions?
    //Trajectory<float_> Traj = Table.getTrajectory(P, 2);

    //std::cout << Traj.getLiftedPositions_x() << std::endl;

    std::vector<float_> Time = seq<float_>(0, 2000, 5000);
    size_t M = 1000000;
    std::vector<float_> EscapeTime(M, 0);
    std::vector<float_> Angle(M, 0);
    for (size_t i = 0; i < M; ++i) {

	Random<float_> StartPos, StartAngle;
	float_ y0 = StartPos.getUniformRandom(Start.P.y, Start.Q.y);
	float_ theta0 = StartAngle.getUniformRandom(0, 2*PI);
	while (cos(theta0) < 0) {
	    theta0 = StartAngle.getUniformRandom(0, 2*PI);
	}
	Angle[i] = theta0;

	Particle<float_> P(Point<float_>(0.000001,y0), Point<float_>(cos(theta0), sin(theta0)), -1, (float_)1e9);
	Particle<float_> P0 = P;

    	//Trajectory<float_> Traj = Table.getTrajectory(P, Time);

	size_t nUpdates = 0;
	while (P.SegmentID != (Lines.size()-1)) {
	    Table.Update(P);
	    nUpdates++;
	}
	EscapeTime[i] = P.TotalTime;
	//std::cout << "Escape time = " << P.TotalTime << std::endl;

    	//Trajectory<float_> Traj = Table.getTrajectory(P0, nUpdates+1);
        //WriteVector("Data/x.dat", Traj.getLiftedPositions_x());
        //WriteVector("Data/y.dat", Traj.getLiftedPositions_y());
    }

    WriteVector("Data/ET.dat", EscapeTime);
    WriteVector("Data/Angle.dat", Angle);

   return 0;
}

