#ifndef _DMWT_H_
#define _DMWT_H_
#include "tetgen.h"
#include "Point3.h"
#include "Vector3.h"
#include "EdgeInfo.h"
#include "TriangleInfo.h"
#include "Configure.h"
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <stdio.h>

#define RIGHT 1  //[a,b]
#define LEFT 0	 //[-infinity,a)&(b,+infinity]
#define triangle(t,i) tris[t*3+i]
#define point(v) Point3(points[v*3],points[v*3+1],points[v*3+2])
#define leftEdgeInd(ei,i) edgeInfoList[ei]->leftEdgeInd[i]
#define leftTris(ei,i) edgeInfoList[ei]->leftTris[i]
#define rightEdgeInd(ei,i) edgeInfoList[ei]->rightEdgeInd[i]
#define rightTris(ei,i) edgeInfoList[ei]->rightTris[i]

class DMWT{
public:
	DMWT();
	DMWT(int ptn, double * pts, double * deGenPts, bool isdegen);
	DMWT(int ptn, double * pts, double * deGenPts, float * norms, bool isdegen);
	~DMWT();
	void preprocess();
	bool start();
	void setWeights(float wtri, float wedge, float wbitri, float wtribd, float wwst);

	// Set triangulation directly (bypasses genTriCandidates/TetGen)
	void setTriangulation(int* triIndices, int numTris);
	
	// Build edge/triangle lists (normally called by preprocess, but needs to be public for setTriangulation)
	void buildList();

	bool reLoadNormalFile(const char* normalfile);
	void statistics();
	void setRound(int r);
	float optimalCost;
	void setDot(bool isdot1);
	void clearTiling();
	void setPointLimit(int limit);

	//------------- for cycle breaking project --------------//
	// outFaces:  flattened triangle vertex positions, size = numTris * 9  (interface 1)
	// outFaces:  flattened triangle vertex indices,  size = numTris * 3   (interface 2)
	// outPoints: curve points, size = numPoints * 3
	// outNorms:  per-point normals, size = numPoints * 3; empty if withNormal == false
	void getResult(std::vector<double>& outFaces, std::vector<double>& outPoints,
	               std::vector<float>& outNorms, bool dosmooth, int subd, int laps);
	void getResult(std::vector<int>& outFaces, std::vector<double>& outPoints,
	               std::vector<float>& outNorms);

	bool EXPSTOP;

protected:
	//-------------variables--------------//
	const char* filename;
	tetgenio in;
	tetgenio out;
	int round;
	int startEdge;
	bool withNormal;
	bool useBiTri;
	bool hasIntersect;//Intersect
	bool hasIntersect2;
	char dot;
	int DMWT_LIMIT;

	int numofpoints;
	int numoftris;
	int numofedges;
	int numofnormals;
	int numoftilingtris;
	int* tris;
	bool owns_tris;  // true if tris was allocated by us (planar path), false if from TetGen
	double* points;
	double* deGenPoints;
	float* normals;
	EdgeInfo ** edgeInfoList;
	TriangleInfo ** triangleInfoList;
	int* tiling;

	float weightTri;
	float weightEdge;
	float weightBiTri;
	float weightTriBd;
	bool useWorstDihedral;

	int ** ehash;
	int ** ehashLeft;
	int ** ehashRight;

	int intsTriInd[2];

	//-------------functions--------------//
	int scanTrianglesOnce();
	char getSide(int v1, int v2, int v3);

	float measureEdge(int v1, int v2);
	float measureTri(int v1, int v2, int v3);
	float measureBiTri(int v1, int v2, int p, int q);
	float measureTriBd(int v1, int v2, int v3, int ni);
	float costTri(float measure);
	float costEdge(float measure);
	float costBiTri(float measure);
	float costTriBd(float measure);

	bool tileSegment(int eind, char side, int ti, float& thiscost, int& thistile);
	void buildTiling(int eind, char side, int ti);
	void genTriCandidates();
	char getSide(int i);

	void initBasics();
	bool triShareEdge(int trii,int trij); 

	// ------------------- for cycle project -----------------//
	bool isDeGen; // degenerated cases: plane
	void saveTilingObj(char* tilefile, const double* finalPoints);
	void saveMeshObj(char* tilefile, int nT, const double* mesh);

	//-------------evaluations--------------//
	float timeReadIn;
	float timePreprocess;
	float timeMWT;
	float timeTotal;
	float timeTetgen;
	float getSize();

	//Not used for now
	int numofcurves;
	bool isOpen;
	int capacity;
	int maxFacePerEdge;

};


#ifndef CURVE_SURFACING_STATIC
#  include "DMWT.cpp"
#  include "DMWT_dot.cpp"
#endif

#endif