#ifndef WRITER_H_INCLUDED
#define WRITER_H_INCLUDED

#include <iostream>
#include <iomanip>
#include <vector>
#include "JumpEvent.h"
#include "LaminarJump.h"

// @brief Wrapper class for writing data (particularly vectors) to a text file
// @app Boolean to determine if data should be appended to file
// @fname String indicating location and name of file (full path)
template <typename T>
class Writer
{

    public:
	Writer() = default;
	//Writer() {

	//    WriteStream.open("default_delme.txt", std::ios::out | std::ios::app);
	//    WriteStream << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;
	//};
	//~Writer(){};
        //Writer(const Writer&) = delete;
        //Writer(Writer&&) = default;

        Writer(const std::string& fname, const bool& app) : FileName(fname), Append(app)
        {
	    if (Append)
		WriteStream.open(FileName.c_str(), std::ios::out | std::ios::app);
	    else 
		WriteStream.open(FileName.c_str(), std::ios::out);
	    WriteStream << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;
        }

	void WriteRowVector(const std::vector<T>& v) {
    	    for (size_t i = 0; i < v.size(); ++i)
    	        WriteStream << v[i] << " ";
    	    WriteStream << std::endl;
	}

	void WriteRowVector(const std::vector<T>& v, const std::vector<size_t>& idx) {
    	    for (size_t i = 0; i < idx.size(); ++i) {
    	        WriteStream << v[idx[i]] << " ";
	    }
    	    WriteStream << std::endl;
	}

	void WriteColVector(const std::vector<T>& v) {
    	    for (size_t i = 0; i < v.size(); ++i)
    	        WriteStream << v[i] << std::endl;
	}

	void WriteColVectors(const std::vector<T>& v1, const std::vector<T>& v2) { 
	    assert(v1.size() == v2.size());
    	    for (size_t i = 0; i < v1.size(); ++i)
    	        WriteStream << v1[i] << " " << v2[i] << std::endl;
	}

	void WriteColVectors(const std::vector<T>& v1, const std::vector<T>& v2, const std::vector<size_t>& v3) {
	    assert(v1.size() == v2.size());
	    assert(v2.size() == v3.size());
    	    for (size_t i = 0; i < v1.size(); ++i)
    	        WriteStream << v1[i] << " " << v2[i] << " " << v3[i] << std::endl;
	}

	//void WriteColVectors(const std::vector<T>& v1, const std::vector<T>& v2, const std::vector<T>& v3) {
	//    assert(v1.size() == v2.size());
	//    assert(v2.size() == v3.size());
    	//    for (size_t i = 0; i < v1.size(); ++i)
    	//        WriteStream << v1[i] << " " << v2[i] << " " << v3[i] << std::endl;
	//}

	void WriteColVectors(const std::vector<T>& v1, const std::vector<T>& v2, const std::vector<T>& v3, const std::vector<size_t>& v4) {
	    assert(v1.size() == v2.size());
	    assert(v2.size() == v3.size());
	    assert(v3.size() == v4.size());
    	    for (size_t i = 0; i < v1.size(); ++i)
    	        WriteStream << v1[i] << " " << v2[i] << " " << v3[i] << " " << v4[i] << std::endl;
	}

	void WriteColVector(const std::vector<std::vector<T>>& v) {

    	    for (size_t i = 0; i < v[0].size(); ++i) { // Loop over outer dimension
		for (size_t j = 0; j < v.size(); ++j) {
		    WriteStream << v[j][i] << " ";
		}
		WriteStream << std::endl;
	    }
	}

	void WriteColVector(const JumpEvent<T>& JE) {
	    WriteStream << JE;
	}

	void Write(const LaminarJump<T>& LJ) {
	    WriteStream << LJ;
	}

	void WriteColElement(const T& v) {
	    WriteStream << v << std::endl;
	}

	void WriteNewLine() {
    	    WriteStream << std::endl;
	}

    private:
	std::ofstream WriteStream;
	std::string FileName;
	bool Append;

};

template <typename T>
void WriteColVector(const std::string &Destination, const JumpEvent<T>& JE) {
    std::ofstream Writer(Destination, std::ios::out | std::ios::app);
    Writer << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;
    Writer << JE;
}

// @brief Writes a vector in one step (there is no point creating an object for this)
template <typename T>
void WriteVector(const std::string &Destination, const std::vector<T> &v)
{
    std::ofstream Writer(Destination, std::ios::out | std::ios::app);
    Writer << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;
    size_t sz = v.size();
    for (size_t i = 0; i < sz; ++i)
        Writer << v[i] << " ";
    Writer << std::endl;
    Writer.close();
}

template <typename T>
void WriteColVector(const std::string &Destination, const std::vector<T> &v)
{
    std::ofstream Writer(Destination, std::ios::out);
    Writer << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;
    size_t sz = v.size();
    for (size_t i = 0; i < sz; ++i)
        Writer << v[i] << std::endl;
    Writer.close();
}

template <typename T>
void WriteColVector(const std::string &Destination, const std::vector<Point<T>> &v)
{
    std::ofstream Writer(Destination, std::ios::out);
    Writer << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;
    size_t sz = v.size();
    for (size_t i = 0; i < sz; ++i)
        Writer << v[i].x << " " << v[i].y << std::endl;
    Writer.close();
}



template <typename T>
void WriteColVector_app(const std::string &Destination, const std::vector<T> &v)
{
    std::ofstream Writer(Destination, std::ios::out | std::ios::app);
    Writer << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;
    for (size_t i = 0; i < v.size(); ++i)
        Writer << v[i] << std::endl;
}

template <typename T>
void WriteParticles(const std::string &Destination, const std::vector<Particle<T>> &P)
{
    std::ofstream Writer(Destination, std::ios::out | std::ios::app);
    Writer << std::setprecision(std::numeric_limits<T>::digits10) << std::showpoint;
    for (size_t i = 0; i < P.size(); ++i) 
	Writer << P[i].Pos.x << " " << P[i].Pos.y << " " << P[i].Vel.x << " " << P[i].Vel.y << std::endl;
    Writer.close();
}


#endif
