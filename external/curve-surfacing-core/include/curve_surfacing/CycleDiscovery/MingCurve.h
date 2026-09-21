#ifndef _MINGCURVE_H_
#define _MINGCURVE_H_

#include "tetgen.h"
#include "Point3.h"
#include "Vector3.h"
#include "Configure.h"
#include "float.h"
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <Eigen/Dense>
#include <igl/fit_plane.h>
#include <igl/triangle/triangulate.h>

using namespace std;

class MingCurve {
public:
	MingCurve(const char* file, int limit, bool hasNorm);
	MingCurve(const double* inCurve,const int inNum, int limit, bool hasNorm);
	MingCurve(const double* inCurve,const float* inNrom, const int inNum, int limit, bool hasNorm);
	~MingCurve();

	int getNumOfPoints();
	double* getPoints();
	double* getDeGenPoints();
	char* getFilename();
	float* getNormal();
	
	// Get 2D triangulation results if available
	bool hasPlanarTriangulation() const { return has_planar_triangulation; }
	const Eigen::MatrixXi& getPlanarTriangulation() const { return planar_triangulation; }

	//----------------Edge protection----------------//
	bool edgeProtect(bool isdmwt);
	void saveCurve(const char* curvefilein, double* pts, int num);

	//-------------evaluations--------------//
	float timeReadIn;
	float timeEdgeProtect;
	void statistics();

	// ------------------- for cycle project -----------------//
	bool isDeGen; // degenerated cases: plane
	bool badInput;

private:
	char* filename;
	int numofpoints;
	double* points;
	double* DeGenPoints;
	float* normals;
	int PT_LIMIT;
	bool EXPSTOP;
	bool withNorm;
	//----------------for edge protection------------//
	bool readOrgCurveFile(const char* file);
	bool loadOrgCurve(const double* inCurve,const int inNum);
	void readOrgNormFile();
	void loadOrgNorm(const float* inNorm);
	int org_n;
	int n_before;
	int n_after;
	float n_ratio;
	std::vector<Point3> tempC;
	std::vector<Point3> tempOrgC;
	std::vector<std::vector<int>> tempAdj;
	std::vector<Vector3> tempNorm;
	std::vector<std::vector<int>> tempAdjNorm;
	void protectCorner();
	double getAngle(int p1, int p2, int p3);
	double getPt2LineDist(int p1, int p2, int p3);
	void splitEdge(int p1, int p2ind, const Point3 & newP);
	void insertMidPointsTetgen();
	std::vector<std::vector<int>> badEdge;
	std::vector<int> newEdge;
	std::vector<int> newAdj;
	std::vector<int> newNorm;
	std::vector<char> newClip;
	bool isProtected();
	void callTetgen();
	
	// ------------------- for cycle project -----------------//
	bool isDeGenCase();
	void perturbPts(double ptb);
	//bool isDeGen; // degenerated cases: plane
	void splitEdge(int p1, int p2ind, const Point3 &newP, const Point3 &newOrgP);
	
	bool getCurveAfterEP();
	bool sameOrientation(const vector<int> & newCurve);
	bool passTetGen();

	// ------------------- 2D planar handling -----------------//
	double computePlanarity();
	bool edgeProtectPlanar2D();
	bool isProtected2D(Eigen::MatrixXd& V2D);  // V2D may be updated to match triangulator output; E is built internally
	void projectTo2D(Eigen::MatrixXd& V2D, Eigen::Vector3d& origin, Eigen::Matrix3d& basis);
	void projectFrom2D(const Eigen::MatrixXd& V2D, const Eigen::Vector3d& origin, const Eigen::Matrix3d& basis);
	
	// Store 2D triangulation results
	Eigen::MatrixXi planar_triangulation;
	bool has_planar_triangulation;

	std::vector<double> radius;
	std::vector<double> orgradius;
	std::vector<std::vector<char>> cliped;
	int perturbNum;

	//not used for now
	int numofcurves;
	int numofnormals;
	bool isOpen;
	int capacity;
};


#ifndef CURVE_SURFACING_STATIC
#  include "MingCurve.cpp"
#endif

#endif