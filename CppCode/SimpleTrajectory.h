#ifndef SIMPLETRAJECTORY_H_INCLUDED
#define SIMPLETRAJECTORY_H_INCLUDED

#include <Geometry.h>

template <typename T>
struct SimpleTrajectory
{
    SimpleTrajectory(const Point<T> &initial_pos, const Point<T>& Diam)
    {
	Diameter = Diam;
	xPos.reserve(1000000);
	Time.reserve(1000000);
	Vel.reserve(1000000);
	cid.reserve(1000000);
        //xPos.push_back(initial_pos.x);
	//Time.push_back(0);
    }

    SimpleTrajectory(const Point<T>& Diam)
    {
	Diameter = Diam;
	//yPos.reserve(10000000);
	xPos.reserve(1000000);
	Time.reserve(1000000);
	Vel.reserve(1000000);
	cid.reserve(1000000);
        //xPos.push_back(initial_pos.x);
	//Time.push_back(0);
    }


    void Push(const Particle<T>& P) {
        xPos.push_back(P.Pos.x + (T)P.CellID.x * Diameter.x);
        //yPos.push_back(P.Pos.y + (T)P.CellID.y * Diameter.y);
	//Time.push_back(P.TotalTime);
	Time.push_back(P.FreePathLength);
	Vel.push_back(P.Vel);
	cid.push_back(P.CellID);
    }

    std::vector<T> xPos;
    //std::vector<T> yPos;
    std::vector<T> Time;
    std::vector<Point<T>> Vel;
    std::vector<Point<long int>> cid;
    Point<T> Diameter;

};


#endif 
