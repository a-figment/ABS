#ifndef BILLIARD_H_INCLUDED
#define BILLIARD_H_INCLUDED

#include <boost/math/constants/constants.hpp>
#include <iomanip>
#include <vector>
#include "Random.h"
#include "Polygon.h"
#include "FloatingPolygon.h"
#include "Circle.h"
#include "ArcSegment.h"
#include "Trajectory.h"
#include "SimpleTrajectory.h"
#include "Beam.h"
#include "Writer.h"

// @brief Is a wrapper of the Obstacle class. Handles state updates for particle<T>s
//

// struct FailState {};

template <typename T>
class Billiard
{

    public:
	Billiard() = default;
	~Billiard(){};

	// Known issues - Even when periodic ID's are set by the obstacle, @Billiard wants to see them explicitly. Fix this by using an obstacle which has meaningful PID parameter
	Billiard(std::vector<std::unique_ptr<Obstacle<T>>> &obstacles)
	{
	    Failed = false;
	    nObstacles = obstacles.size();
	    Obstacles = std::move(obstacles);
	    Diameter = Point<T>(Obstacles[0]->getCellEdge().first.y - Obstacles[0]->getCellEdge().first.x, Obstacles[0]->getCellEdge().second.x - Obstacles[0]->getCellEdge().second.y);
	    Diameter.y = 0.;
	    std::cout << "WARNING: Vertical diameter is being forced in Billiard.h" << std::endl;
	    //std::cout << Obstacles[0]->getCellEdge().second.x - Obstacles[0]->getCellEdge().second.y << std::endl;
	    // The zeroth obstacle MUST contain information about periodic boundaries
	    LRpid = Obstacles[0]->getPeriodicIDs()[0]; 
	    UDpid = Obstacles[0]->getPeriodicIDs()[1]; 


	}

	void Update(Particle<T> &P);
	void UpdateCollisionState(Particle<T> &P);
	void UpdatePeriodicState(Particle<T> &P);

	Particle<T> getInitialState(const int &OID, const int &SID, const T &mt) const;
	Particle<T> getIC_FixedPosition(const int &OID, const int &SID, const T &mt, const Point<T> &Position) const;
	Particle<T> getIC_FixedVelocity(const int &OID, const int &SID, const T &mt, const T &Angle) const;
	Particle<T> getIC_FixedLine(const int &OID, const int &SID, const T &mt, const Point<T>& A, const Point<T>& B) const;

	Particle<T> getIC_Rectangle(const int &OID, const int &SID, const T &mt, const Point<T>& A, const Point<T>& B) const;
	Particle<T> getIC_ExpandingBall(const int &OID, const int &SID, const T &mt, const Point<T>& center, const T& radius);

	// @brief
	Trajectory<T> getTrajectory(Particle<T> &P, const std::vector<T> &Time);
	Trajectory<T> getTrajectory(Particle<T> &P, const size_t &nCollisions);
	Trajectory<T> getTrajectory(Particle<T> &P, const size_t &nCollisions, const T& mt, const T& Mt);
	SimpleTrajectory<T> getSimpleTrajectory(Particle<T> &P);
	void WriteLongTimeSeries(Particle<T> &P, const size_t &nCollisions, const std::string& DataDir);

	// @brief Finite trajectories
	std::pair<Trajectory<T>, std::vector<int>> getFiniteTrajectory(Particle<T> &P, const std::vector<T> &Time, const std::vector<T> &L); // Continuous and finite
	Trajectory<T> getFiniteTrajectory(Particle<T> &P, const size_t &nCollisions, const T &L);
	Trajectory<T> getFiniteTrajectory(Particle<T> &P, const T &L);

	// Unforunately FiniteTrajectories use a lot of memory so we need break it up.
	std::pair<std::vector<T>, std::vector<int>> getEscapeTimes(Particle<T> &P, const std::vector<T> &L);

	// @brief Gets discrete and continuous particle<T> state
	Trajectory<T> getTrajectorySample(Particle<T> &P, const std::vector<T> &Time);

	Beam<T> getInitialBeam();
	Beam<T> getInitialBeamWithPB();
	// std::pair<Obstacle<T>, size_t> getExternal();
	// std::unique_ptr<Obstacle<T>> getExternal();
	bool Error(const Particle<T> &P)
	{
	    Failed = !P.collisionFound;
	    return Failed;
	}

	void ErrorReport(const Particle<T> &P, const size_t &index)
	{
	    std::cout << ">>>>>>>>>>>>>>>>> ERROR REPORT: COLLISION NOT FOUND <<<<<<<<<<<<<<<<<" << std::endl;
	    std::cout << "Particle index    = " << index << "," << std::endl;
	    std::cout << "Number of updates = " << P.nUpdates << "." << std::endl;
	    std::cout << " >>>>> STATE AT ERROR: " << std::endl;
	    std::cout << P << std::endl
		<< std::endl;

	    // std::cout << " >>>>> STATE PRIOR TO ERROR: " << std::endl;
	    // Particle<T> Q(P.Pos, Point<T>(-P.Vel.x, -P.Vel.y), P.SegmentID, P.MaxTime);
	    // UpdateCollisionState(Q);
	    // std::cout << Q << std::endl;
	}
	void UpdatePBCs(Particle<T> &P);

	//T getMFPT() const { return this->MFPT; }
	void setMFPT();

	// Point<T> getPointInside(const Billiard &D) const;
	std::vector<std::unique_ptr<Obstacle<T>>> Obstacles;
	bool Failed; // A particular state could not be updated

	// @brief Returns points (s, phi) = (arc length (CCW), angles) [CURRENTLY WORKS ONLY FOR SINGLE OBSTACLE BILLIARDS]
	std::vector<Point<T>> getBirkhoffCoords(const Trajectory<T> &Traj) const;

	Point<T> getDiameter() const { return Diameter; }

    private:
	bool onPBC(const Particle<T> &P);
	size_t nObstacles;
	T MFP;
	Point<T> Diameter;
	std::pair<int, int> LRpid;
	std::pair<int, int> UDpid;
};

// @brief Interpolates along the arc length of the billiard path in steps of dt = Time[i] - Time[i-1]
    template <typename T>
Trajectory<T> Billiard<T>::getTrajectory(Particle<T> &P, const std::vector<T> &Time)
{
    Trajectory<T> Traj(Time.size(), P, Diameter);
    size_t collCount = 1, timeCount = 1;
    while (timeCount < Time.size() && P.TotalTime <= Time.back())
    {                      // Get collisions until path length reaches T
	Particle<T> tP = P; // Particle<T> in continuous time
	Update(P);          // Find the next collision
	collCount++;
	if (Error(P))
	{
	    ErrorReport(tP, collCount);
	    break;
	}

	// Interpolate along collision trajectory
	while (timeCount < Time.size() && P.TotalTime > Time[timeCount])
	{
	    T tau = Time[timeCount] - tP.TotalTime;
	    tP.UpdateState(tau, -1, -1); // Continuous trajectory assumed to never lie on an obstacle
	    Traj.Push(tP);
	    //std::cout << tP << std::endl;
	    timeCount++;
	}
    }
    Traj.nCollisions = collCount;
    return Traj;
}

// @brief Trajectory in terms of collisions
    template <typename T>
Trajectory<T> Billiard<T>::getTrajectory(Particle<T> &P, const size_t &nCollisions)
{
    Trajectory<T> Traj(nCollisions, P, Diameter);
    for (size_t i = 0; i < nCollisions - 1; ++i)
    {
	Particle<T> P0 = P;
	Update(P); // Find the next collision

	if (Error(P))
	{
	    std::cout << "STATE PRIOR TO ERROR: " << std::endl;
	    std::cout << P0 << std::endl;
	    P0.WriteState("Data/ErrorState.dat");
	    break;
	}
	Traj.Push(P);
    }
    return Traj;
}

// @brief Trajectory in terms of collisions up to a maximum total time @P.MaxTime
    template <typename T>
SimpleTrajectory<T> Billiard<T>::getSimpleTrajectory(Particle<T> &P)
{
    //SimpleTrajectory<T> STraj(P.Pos, Diameter);
    //T t = STraj.Time.back(); // For continuing a trajectory
    SimpleTrajectory<T> STraj(Diameter);
    T t = (T)0;

    while(t < P.MaxTime)
    {
	STraj.Push(P);
	Update(P); // Find the next collision

	t = P.TotalTime;
    }
    return STraj;
}


// @brief Writes (t,x) data. Useful when @nCollisions is very large
    template <typename T>
void Billiard<T>::WriteLongTimeSeries(Particle<T> &P, const size_t &nCollisions, const std::string& DataDir)
{
    Writer<T> _X(DataDir + "x.dat", true);
    Writer<T> _t(DataDir + "t.dat", true);

    Point<T> prevPos(P.Pos.x, P.Pos.y);
    T prevLen = 0;

    _X.WriteColElement(prevPos.x);
    _t.WriteColElement(prevLen);

    for (size_t i = 1; i < nCollisions; ++i)
    {
	Update(P); // Find the next collision

	Point<T> currPos = P.Pos + Point<T>((T)P.CellID.x * Diameter.x, (T)P.CellID.y * Diameter.y); 
	T currLen = Distance2D<T>(prevPos, currPos);	

	_X.WriteColElement(currPos.x);
	_t.WriteColElement(currLen);

	prevPos = currPos;
	prevLen = currLen;
    }
}


// @brief Gets trajectory by collision such that the length is between min time @mt and max time @Mt
// A large @nCollisions is recommended since
    template <typename T>
Trajectory<T> Billiard<T>::getTrajectory(Particle<T> &P, const size_t &nCollisions, const T& mt, const T& Mt)
{
    Trajectory<T> Traj(nCollisions, Diameter);
    while (P.TotalTime < mt)
    {
	Particle<T> P0 = P;
	Update(P); // Find the next collision

	if (Error(P))
	{
	    std::cout << "STATE PRIOR TO ERROR: " << std::endl;
	    std::cout << P0 << std::endl;
	    P0.WriteState("Data/ErrorState.dat");
	    break;
	}
    }

    //std::cout << "Start record time = " << P.TotalTime << std::endl;

    // Start recording states onces TotalTime > min time
    while (P.TotalTime < Mt)
    {
	Particle<T> P0 = P;
	Update(P); // Find the next collision

	if (Error(P))
	{
	    std::cout << "STATE PRIOR TO ERROR: " << std::endl;
	    std::cout << P0 << std::endl;
	    P0.WriteState("Data/ErrorState.dat");
	    break;
	}
	Traj.Push_(P);
    }

    //std::cout << "End record time = " << P.TotalTime << std::endl;

    return Traj;
}


// @brief Collision trajectory up to a finite distance L.
// @return A trajectory of nCollisions which did not escape or one of < nCollisions which did.
    template <typename T>
Trajectory<T> Billiard<T>::getFiniteTrajectory(Particle<T> &P, const size_t &nCollisions, const T &L)
{
    Trajectory<T> Traj(nCollisions, P, Diameter);

    for (size_t i = 0; i < nCollisions - 1; ++i)
    {
	// Particle<T> oldP = P;
	Update(P); // Find the next collision
	if (Error(P) || P.TotalTime > P.MaxTime)
	    break;
	Traj.Push(P);

	T X = P.Pos.x + static_cast<T>(P.CellID.x);
	if (fabs(X) > L) // Escape
	{
	    Traj.Pop(Traj.getNumElements()); // Returns trajectory with final state occurring immediately after an escape
	    break;
	}
    }
    return Traj;
}

// @brief Collision trajectory up to a finite distance L OR total time P.MaxTime
// @return
    template <typename T>
Trajectory<T> Billiard<T>::getFiniteTrajectory(Particle<T> &P, const T &L)
{
    Trajectory<T> Traj; //(static_cast<size_t>(P.MaxTime), P);
    Traj.Reserve(static_cast<size_t>(P.MaxTime), Diameter);
    Traj.Push_(P);

    T X = P.Pos.x + static_cast<T>(P.CellID.x);
    while (fabs(X) < L && P.TotalTime <= P.MaxTime)
    {
	Update(P); // Find the next collision
	X = P.Pos.x + static_cast<T>(P.CellID.x);
	if (Error(P))
	    break;
	Traj.Push_(P);
	// Traj.Push(P);
    }
    return Traj;
}

// @brief Finite continuous trajectory
    template <typename T>
std::pair<Trajectory<T>, std::vector<int>> Billiard<T>::getFiniteTrajectory(Particle<T> &P, const std::vector<T> &Time, const std::vector<T> &L)
{
    Trajectory<T> Traj(Time.size(), P, Diameter);
    std::vector<int> EscapeIndices(L.size(), -1);

    size_t collCount = 1, timeCount = 1;
    for (size_t i = 0; i < L.size(); ++i)
    {
	while (timeCount < Time.size() && P.TotalTime <= Time.back())
	{                      // Get collisions until path length reaches T
	    Particle<T> tP = P; // Particle<T> in continuous time
	    Update(P);          // Find the next collision
	    collCount++;
	    if (Error(P))
		break;

	    // Interpolate along collision trajectory
	    while (timeCount < Time.size() && P.TotalTime > Time[timeCount])
	    {
		T tau = Time[timeCount] - tP.TotalTime;
		tP.UpdateState(tau, -1, -1); // Continuous trajectory assumed to never lie on an obstacle
		Traj.Push(tP);
		timeCount++;
	    }

	    T X = P.Pos.x + (T)P.CellID.x;
	    if (fabs(X) > L[i])
	    {
		EscapeIndices[i] = (int)timeCount;
		break;
	    }
	}
    }

    // if ((size_t)EscapeIndices.back() < Time.size())
    //    Traj.Pop(Time.size() - (int)EscapeIndices.back() + 1);

    return std::pair<Trajectory<T>, std::vector<int>>(Traj, EscapeIndices);
}

// @brief Gets the time of escape when x(t) = L. Less memory intensive if large times are needed.
template <typename T>
std::pair<std::vector<T>, std::vector<int>>  Billiard<T>::getEscapeTimes(Particle<T> &P, const std::vector<T> &L)
{
    std::vector<T> FPT(L.size(), INFINITY);
    std::vector<int> nCollisions(L.size(), 0);
    T X = P.Pos.x + static_cast<T>(P.CellID.x);

    for (size_t i = 0; i < L.size(); ++i)
    {
	size_t nCollisionsAtEscape = 0;
	Particle<T> PriorToEscape = P;
	while (fabs(X) < L[i] && P.TotalTime <= P.MaxTime)
	{
	    PriorToEscape = P;
	    Update(P); // Find the next collision
	    X = P.Pos.x + static_cast<T>(P.CellID.x);
	    nCollisionsAtEscape++;

	    if (Error(P))
		break;
	}
	if (fabs(X) > L[i])
	{ // An escape occurred
	    T d = getDistanceToEscape<T>(PriorToEscape, L[i]);
	    FPT[i] = PriorToEscape.TotalTime + d;

	    if (X > 0)
		nCollisions[i] = P.nUpdates;
	    else
		nCollisions[i] = -P.nUpdates;
	    P = PriorToEscape; // For infinite horizons
	}
	else if (P.TotalTime > P.MaxTime)
	{
	    return std::make_pair(FPT, nCollisions); //
	    std::cout << "Exceeded max time for L = " << L[i] << std::endl;
	}
    }
    return std::make_pair(FPT, nCollisions);
}

// @brief Gets a trajectory of length @Time.back() with both collisions and interpolated points
    template <typename T>
Trajectory<T> Billiard<T>::getTrajectorySample(Particle<T> &P, const std::vector<T> &Time)
{
    Trajectory<T> Traj;
    Traj.Reserve(2 * Time.size(), Diameter);
    Traj.Push_(P);

    size_t colCount = 1, timeCount = 1;
    while (P.TotalTime <= Time.back() && timeCount < Time.size())
    {                      // Get collisions until path length reaches T
	Particle<T> tP = P; // Particle<T> in continuous time
	Update(P);          // Find the next collision
	colCount++;
	if (Error(P) || P.TotalTime > P.MaxTime)
	    break;

	// Interpolate along collision trajectory
	while (P.TotalTime > Time[timeCount] && timeCount < Time.size())
	{
	    T tau = Time[timeCount] - tP.TotalTime;
	    tP.UpdateState(tau, -1, -1); // Continuous trajectory assumed to never lie on an obstacle
	    Traj.Push_(tP);
	    timeCount++;
	}
	Traj.Push_(P);
    }
    return Traj;
}

// @brief Returns points (s, phi) = (arc length (CCW), angles)
// CURRENTLY WORKS ONLY FOR SINGLE OBSTACLE BILLIARDS
template <typename T>
std::vector<Point<T>> Billiard<T>::getBirkhoffCoords(const Trajectory<T> &Traj) const
{
    std::vector<T> SIDs = Traj.getSegmentIDs();
    std::vector<T> OIDs = Traj.getObstacleIDs();
    std::vector<Point<T>> Positions = Traj.getPositions();
    std::vector<Point<T>> Velocities = Traj.getVelocities();
    std::vector<Point<T>> Vertices = Traj.getVelocities();

    std::vector<Point<T>> BCs(SIDs.size(), Point<T>(-1, -1)); // Birkhoff coordinates
    for (size_t i = 1; i < BCs.size(); ++i)
    { // First position may be on the boundary

	Point<T> Vertex = Obstacles[OIDs[i]]->getVertices()[SIDs[i]];
	Point<T> Normal = Obstacles[OIDs[i]]->getInwardNormal(Positions[i], SIDs[i]);

	T s = Distance2D<T>(Positions[i], Vertex) + Obstacles[OIDs[i]]->getCumulativePerimiter()[SIDs[i]];
	T phi = acos(Dot<T>(Normal, Velocities[i]));

	BCs[i] = Point<T>(s, phi);
    }
    return BCs;
}

    template <typename T>
void Billiard<T>::UpdateCollisionState(Particle<T> &P)
{
    P.collisionFound = false; // Remains false if P is not updated

    if (Obstacles.size() == 1) // Collide finds the closest collision
    {
	Obstacles[0]->Collide(P, static_cast<int>(0));
    }
    else
    {
	T CollisionTime = INFINITY;
	Particle<T> NewState;
	for (size_t i = 0; i < nObstacles; ++i)
	{
	    Particle<T> OldState = P;
	    Obstacles[i]->Collide(OldState, static_cast<int>(i));

	    if (OldState.collisionFound)
	    {
		if (OldState.FreePathLength < CollisionTime)
		{
		    NewState = OldState;
		    CollisionTime = OldState.FreePathLength;
		}
	    }
	}
	P = NewState;
    }
}

// @brief Finds the next collision which is not on a periodic boundary
template <typename T>
void Billiard<T>::Update(Particle<T> &P)
{
    UpdateCollisionState(P);
    if (!P.collisionFound)
    {
	std::cout << "ERROR: Collision not found. Please use Billiard::ErrorReport() for more info." << std::endl;
	return;
    }
    //std::cout << P << std::endl;
    UpdatePeriodicState(P);
}

// @brief Finds the next collision which is not on a periodic boundary
    template <typename T>
void Billiard<T>::UpdatePeriodicState(Particle<T> &P)
{
    UpdatePBCs(P);


    while (onPBC(P)) // Handle infinite horizon
    {
	//std::cout << P << std::endl;
	//std::cout << "made it" << std::endl;
	UpdateCollisionState(P);
	P.nUpdates--;
	UpdatePBCs(P);
	if (P.TotalTime > P.MaxTime) // Infinite horizon billiards may stay in this loop for longer than the total time
	    break;
    }
}

template <typename T>
bool Billiard<T>::onPBC(const Particle<T> &P)
{
    // Since the first obstacle is always the exterior, it can be the only one to have PBCs
    //std::pair<int, int> LRpid = Obstacles[0]->getPeriodicIDs()[0]; // LEFT/RIGHT
    //std::pair<int, int> UDpid = Obstacles[0]->getPeriodicIDs()[1]; // UP/DOWN

    return ((P.ObstacleID == 0) && ((P.SegmentID == LRpid.first) || (P.SegmentID == LRpid.second) || (P.SegmentID == UDpid.first) || (P.SegmentID == UDpid.second)));
}

template <typename T>
void Billiard<T>::UpdatePBCs(Particle<T> &P) // Only handles translations in horizontal/vertical directions
{
    //constexpr std::pair<int, int> LRpid = Obstacles[0]->getPeriodicIDs()[0]; 
    //constexpr std::pair<int, int> UDpid = Obstacles[0]->getPeriodicIDs()[1]; 

    if (P.ObstacleID != 0)
    {
	P.Vel = ReflectSpecular<T>(Obstacles[P.ObstacleID]->getInwardNormal(P), P.Vel);
	//std::cout << P.Vel << std::endl;
	//std::cout << Obstacles[P.ObstacleID]->getInwardNormal(P) << std::endl;

	return;
    }
    //std::cout << Obstacles[0]->getCellEdge().first.x << std::endl;


    if (P.SegmentID == LRpid.first) //
    {
	P.SegmentID = LRpid.second;
	P.CellID.x += 1;
	P.Pos.x = Obstacles[0]->getCellEdge().first.x;
    }
    else if (P.SegmentID == LRpid.second)
    {
	P.SegmentID = LRpid.first;
	P.CellID.x -= 1;
	P.Pos.x = Obstacles[0]->getCellEdge().first.y;
    }
    else if (P.SegmentID == UDpid.first)
    {
	P.SegmentID = UDpid.second;
	P.CellID.y -= 1;
	P.Pos.y = Obstacles[0]->getCellEdge().second.x;
    }
    else if (P.SegmentID == UDpid.second)
    {
	P.SegmentID = UDpid.first;
	P.CellID.y += 1;
	P.Pos.y = Obstacles[0]->getCellEdge().second.y;
    }
    else
    {
	P.Vel = ReflectSpecular<T>(Obstacles[P.ObstacleID]->getInwardNormal(P), P.Vel);
    }

    //switch (P.SegmentID) {
    //    case LRpid.first:
    //        P.SegmentID = LRpid.second;
    //        P.CellID.x += 1;
    //  	    P.Pos.x = Obstacles[0]->getCellEdge().first.x;
    //  	    break;
    //  	case LRpid.second:
    //        P.SegmentID = LRpid.first;
    //  	    P.CellID.x -= 1;
    //  	    P.Pos.x = Obstacles[0]->getCellEdge().first.y;
    //        break; 
    //  	case UDpid.first:
    //        P.SegmentID = UDpid.second;
    //  	    P.CellID.y -= 1;
    //  	    P.Pos.y = Obstacles[0]->getCellEdge().second.x;
    //        break;
    //  	case UDpid.second:
    //        P.SegmentID = UDpid.first;
    //  	    P.CellID.y += 1;
    //  	    P.Pos.y = Obstacles[0]->getCellEdge().second.y;
    //  	    break;
    //  	default: 
    //        P.Vel = ReflectSpecular<T>(Obstacles[P.ObstacleID]->getInwardNormal(P), P.Vel);
    //}
}

// ============ GETTERS ============

// @brief
// @id If id == -1, particle<T> position should be floating on the interior of billiards
//     If id != -1, particle<T> may be
template <typename T>
Particle<T> Billiard<T>::getInitialState(const int &OID, const int &SID, const T &mt) const
{

    //TODO - Fix this for disk channel

    // Obstacle[0] is always the exterior obstacle
    size_t counter = 0;
    bool inside = false;
    Point<T> InitialPosition(0, 0);
    while (!inside) // Find a point inside Obstacle[0] but not in Obstacle[i], i > 0
    {
	InitialPosition = Obstacles[0]->getPointInside();

	if (Obstacles.size() > 1)
	{
	    for (size_t i = 1; i < Obstacles.size(); ++i)
	    {
		if (Obstacles[i]->isPointInside(InitialPosition))
		    break;
		else
		    inside = true;
	    }
	}
	else
	{
	    inside = true;
	}
	counter++;
	if (counter > 1000000)
	    std::cout << "Struggling to find initial condition! " << counter << " attempts made." << std::endl;
    }

    // Uniform random angles
    Random<T> URAngle;
    T Theta = URAngle.getUniformRandom(static_cast<T>(0), static_cast<T>(2) * boost::math::constants::pi<T>());
    Point<T> InitialVelocity(cos(Theta), sin(Theta));

    Particle<T> p(InitialPosition, InitialVelocity, OID, SID, mt);
    //p.Diameter = Point<T>(Obstacles[0]->getCellEdge().first.y - Obstacles[0]->getCellEdge().first.x, Obstacles[0]->getCellEdge().second.x - Obstacles[0]->getCellEdge().second.y);

    return p;
}

// Special case...
//template <typename T>
//Particle<T> Billiard<T>::getInitialState_RoundedChannel(const int &OID, const int &SID, const T &mt) const
//{
//
//    Random<T> rx;
//    Random<T> ry;
//
//    Point<T> testPoint;
//    bool isInside = false;
//    while (!isInside) {
//	testPoint = Point<T>(rx.getUniformRandom(Obstacles[0]->BoundingBox.first.x, Obstacles[0]->BoundingBox.second.x), // Base polygon
//			     ry.getUniformRandom(Obstacles[1]->BoundingBox.first.y, this->BoundingBox.second.y));
//
//	
//	
//	isInside = isPointInside(testPoint);
//    }
//
//
//    Random<T> UR;
//    T Theta = URAngle.getUniformRandom(static_cast<T>(0), static_cast<T>(2) * boost::math::constants::pi<T>());
//
//    // Obstacle[0] is always the exterior obstacle
//    size_t counter = 0;
//    bool inside = false;
//    Point<T> InitialPosition(0, 0);
//    while (!inside) // Find a point inside Obstacle[0] but not in Obstacle[i], i > 0
//    {
//	InitialPosition = Obstacles[0]->getPointInside();
//
//	if (Obstacles.size() > 1)
//	{
//	    for (size_t i = 1; i < Obstacles.size(); ++i)
//	    {
//		if (Obstacles[i]->isPointInside(InitialPosition))
//		    break;
//		else
//		    inside = true;
//	    }
//	}
//	else
//	{
//	    inside = true;
//	}
//	counter++;
//	if (counter > 1000000)
//	    std::cout << "Struggling to find initial condition! " << counter << " attempts made." << std::endl;
//    }
//
//    // Uniform random angles
//    Random<T> URAngle;
//    T Theta = URAngle.getUniformRandom(static_cast<T>(0), static_cast<T>(2) * boost::math::constants::pi<T>());
//    Point<T> InitialVelocity(cos(Theta), sin(Theta));
//
//    Particle<T> p(InitialPosition, InitialVelocity, OID, SID, mt);
//    //p.Diameter = Point<T>(Obstacles[0]->getCellEdge().first.y - Obstacles[0]->getCellEdge().first.x, Obstacles[0]->getCellEdge().second.x - Obstacles[0]->getCellEdge().second.y);
//
//    return p;
//}


// @brief Distributes positions uniformly along the line A-B with uniform velocity
template <typename T>
Particle<T> Billiard<T>::getIC_FixedLine(const int &OID, const int &SID, const T &mt, const Point<T>& A, const Point<T>& B) const
{
    // Uniform random angles
    Random<T> URAngle;
    T Theta = URAngle.getUniformRandom(static_cast<T>(0), static_cast<T>(2) * boost::math::constants::pi<T>());
    Point<T> InitialVelocity(cos(Theta), sin(Theta));

    // Uniform points on the line segment
    Random<T> Param;
    T t = Param.getUniformRandom((T)0, (T)1);
    T x = (1-t)*A.x + t*B.x;
    T y = (1-t)*A.y + t*B.y;
    Point<T> Position(x,y);

    Particle<T> p(Position, InitialVelocity, OID, SID, mt);

    return p;
}


// @brief Distributes positions uniformly in a rectangular region with bounding coordinates (@A, @B)
template <typename T>
Particle<T> Billiard<T>::getIC_Rectangle(const int &OID, const int &SID, const T &mt, const Point<T>& A, const Point<T>& B) const
{
    // Uniform random angles
    Random<T> URAngle;
    T Theta = URAngle.getUniformRandom(static_cast<T>(0), static_cast<T>(2) * boost::math::constants::pi<T>());
    Point<T> InitialVelocity(cos(Theta), sin(Theta));

    // Uniform points on the line segment
    Random<T> Px, Py;

    T X = Px.getUniformRandom(A.x, B.x);
    T Y = Py.getUniformRandom(A.y, B.y);

    Point<T> Position(X,Y);
    Particle<T> p(Position, InitialVelocity, OID, SID, mt);

    return p;
}



// @brief Distribute points uniformly on a circle with @center and @radius. Velocity of particle is taken to be the outward normal of the ball
    template <typename T>
Particle<T> Billiard<T>::getIC_ExpandingBall(const int &OID, const int &SID, const T &mt, const Point<T>& center, const T& radius)
{
    // Uniform random angles
    Random<T> URAngle;
    T Theta = URAngle.getUniformRandom(static_cast<T>(0), static_cast<T>(2) * boost::math::constants::pi<T>());
    Point<T> InitialVelocity(cos(Theta), sin(Theta));

    // Uniform points on the line segment
    Random<T> Radius((T)0, (T)radius);
    T r = Radius.getUniformRandom();
    //T phi = URAngle.getUniformRandom();

    Point<T> Position(center.x + sqrt(r) * cos(Theta), center.y + sqrt(r) * sin(Theta));

    Particle<T> p(Position, InitialVelocity, OID, SID, mt);

    return p;
}

// @brief G
//template <typename T>
//Particle<T> Billiard<T>::getIC_MovingBall(const int &OID, const int &SID, const T &mt, Point<T>& center, const T& radius)
//{
//   // Uniform random angles
//   Random<T> URAngle(static_cast<T>(0), static_cast<T>(2) * boost::math::constants::pi<T>());
//   T Theta = URAngle.getUniformRandom();
//   Point<T> InitialVelocity(cos(Theta), sin(Theta));
//  
//   // Uniform points on the line segment
//   Random<T> Param((T)0, (T)1);
//   T t = Param.getUniformRandom();
//   T x = (1-t)*A.x + t*B.x;
//   T y = (1-t)*A.y + t*B.y;
//   Point<T> Position(x,y);
//
//   Particle<T> p(Position, InitialVelocity, OID, SID, mt);
//
//   return p;
//}


template <typename T>
Particle<T> Billiard<T>::getIC_FixedPosition(const int &OID, const int &SID, const T &mt, const Point<T>& Position) const
{
    // Uniform random angles
    Random<T> URAngle;
    T Theta = URAngle.getUniformRandom(static_cast<T>(0), static_cast<T>(2) * boost::math::constants::pi<T>());
    Point<T> InitialVelocity(cos(Theta), sin(Theta));

    Particle<T> p(Position, InitialVelocity, OID, SID, mt);
    //p.Diameter = Point<T>(Obstacles[0]->getCellEdge().first.y - Obstacles[0]->getCellEdge().first.x, Obstacles[0]->getCellEdge().second.x - Obstacles[0]->getCellEdge().second.y);

    return p;
}

template <typename T>
Particle<T> Billiard<T>::getIC_FixedVelocity(const int &OID, const int &SID, const T &mt, const T& Angle) const
{

    // Obstacle[0] is always the exterior obstacle
    bool inside = false;
    Point<T> InitialPosition(0, 0);
    while (!inside) // Find a point inside Obstacle[0] but not in Obstacle[i], i > 0
    {
	InitialPosition = Obstacles[0]->getPointInside();

	if (Obstacles.size() > 1)
	{
	    for (size_t i = 1; i < Obstacles.size(); ++i)
	    {
		if (Obstacles[i]->isPointInside(InitialPosition))
		    break;
		else
		    inside = true;
	    }
	}
	else
	{
	    inside = true;
	}
    }

    Point<T> InitialVelocity(cos(Angle), sin(Angle));

    Particle<T> p(InitialPosition, InitialVelocity, OID, SID, mt);
    //p.Diameter = Point<T>(Obstacles[0]->getCellEdge().first.y - Obstacles[0]->getCellEdge().first.x, Obstacles[0]->getCellEdge().second.x - Obstacles[0]->getCellEdge().second.y);

    return p;
}

// @brief The beam information must be fixed before each 'collection' of particles is considered. We generate a random point
//        and direction in the billiard and cast a ray. The points of intersection determine all positions on the beam. A new angle is generated which
//        lies between the tangent line defined at either one of the two end points.
//        ASSUMES NO INTERIOR OBSTACLES.
// @SID The beam by definition connects two segments of the billiard. In principle it is possible that SEG.x = SEG.y.
// @OID The beam may connect more than one obstacle
// @n The nth particle along the beam
    template <typename T>
Beam<T> Billiard<T>::getInitialBeam() 
{
    // Here we take a valid point inside the billiard as normal.
    const T PI = boost::math::constants::pi<T>();
    bool inside = false;
    Point<T> InitialPosition(0, 0);
    while (!inside) {
	InitialPosition = Obstacles[0]->getPointInside();
	if (Obstacles.size() > 1)
	{
	    for (size_t i = 1; i < Obstacles.size(); ++i)
	    {
		if (Obstacles[i]->isPointInside(InitialPosition))
		    break;
		else
		    inside = true;
	    }
	}
	else
	{
	    inside = true;
	}
    }


    // Generate ray which doesn't pass through the periodic boundary -- this is a conceptual difference more than technical issue. We want to maintain the typical transport propbel of uniform positions/velicites
    Random<T> URAngle;
    Particle<T> RayParticleA;
    Particle<T> RayParticleB;
    T RayAngle;
    bool onPeriodic = true;
    while (onPeriodic) {
	// Uniform random angle
	RayAngle = URAngle.getUniformRandom((T)0, (T)2 * PI);

	// The line passing through InitialPosition at angle RayAngle will intersect the billiard table. 
	RayParticleA = Particle<T>(InitialPosition, Point<T>(cos(RayAngle), sin(RayAngle)), 0, -1, (T)0);
	RayParticleB = Particle<T>(InitialPosition, Point<T>(-cos(RayAngle), -sin(RayAngle)), 0, -1, (T)0);
	UpdateCollisionState(RayParticleA);
	UpdateCollisionState(RayParticleB);
	onPeriodic = onPBC(RayParticleA) || onPBC(RayParticleB); // Not on periodic, so exit
    }
    Point<T> A = RayParticleA.Pos;
    Point<T> B = RayParticleB.Pos;

    Point<int> SIDs(RayParticleA.SegmentID, RayParticleB.SegmentID);

    T DepartureAngle = URAngle.getUniformRandom((T)0, (T)2*PI);
    Beam<T> beam(A, B, RayAngle, DepartureAngle);

    return beam;
}

    template <typename T>
Beam<T> Billiard<T>::getInitialBeamWithPB() 
{
    // Here we take a valid point inside the billiard as normal.
    const T PI = boost::math::constants::pi<T>();
    bool inside = false;
    Point<T> InitialPosition(0, 0);
    while (!inside) {
	InitialPosition = Obstacles[0]->getPointInside();
	if (Obstacles.size() > 1)
	{
	    for (size_t i = 1; i < Obstacles.size(); ++i)
	    {
		if (Obstacles[i]->isPointInside(InitialPosition))
		    break;
		else
		    inside = true;
	    }
	}
	else
	{
	    inside = true;
	}
    }


    // Generate ray which might pass through periodic boundary
    Random<T> URAngle;
    Particle<T> RayParticleA;
    Particle<T> RayParticleB;
    T RayAngle;
    // Uniform random angle
    RayAngle = URAngle.getUniformRandom((T)0, (T)2 * PI);

    // The line passing through InitialPosition at angle RayAngle will intersect the billiard table. 
    RayParticleA = Particle<T>(InitialPosition, Point<T>(cos(RayAngle), sin(RayAngle)), 0, -1, (T)0);
    RayParticleB = Particle<T>(InitialPosition, Point<T>(-cos(RayAngle), -sin(RayAngle)), 0, -1, (T)0);
    UpdateCollisionState(RayParticleA);
    UpdateCollisionState(RayParticleB);
    Point<T> A = RayParticleA.Pos;
    Point<T> B = RayParticleB.Pos;

    Point<int> SIDs(RayParticleA.SegmentID, RayParticleB.SegmentID);

    T DepartureAngle = URAngle.getUniformRandom((T)0, (T)2*PI);
    Beam<T> beam(A, B, RayAngle, DepartureAngle);

    return beam;
}
// @brief Accepts a previously generated beam datum and returns a particle not previously considered
//template <typename T>
//Particle<T> Billiard<T>::getParticleOnBeam() const
//{
//
//   Particle<T> p(InitialPosition, InitialVelocity, OID, SID, mt);
//   p.Diameter = Point<T>(Obstacles[0]->getCellEdge().first.y - Obstacles[0]->getCellEdge().first.x, Obstacles[0]->getCellEdge().second.x - Obstacles[0]->getCellEdge().second.y);
//
//   return p;
//}



// ============ SETTERS ============

// @TODO: Currently incorrect for billiards which may have internal scatterers.
// Obstacles should have an 'internal' or 'external' flag
    template <typename T>
void Billiard<T>::setMFPT()
{
    T TotalArea = static_cast<T>(0);
    T TotalPerimeter = static_cast<T>(0);
    for (size_t i = 0; i < nObstacles; ++i)
    {
	TotalArea += Obstacles[i]->getArea();
	TotalPerimeter += Obstacles[i]->getPerimeter();
    }
    MFP = boost::math::constants::pi<T>() * (TotalArea / TotalPerimeter);
}

#endif //BILLIARD_H_INCLUDED
