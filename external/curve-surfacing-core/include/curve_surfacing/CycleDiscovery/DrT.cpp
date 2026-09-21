#pragma once
#include "../curve_surfacing_inline.h"

// ── implementation includes ──
#include "DrT.h"
#include "../cs_logging.h"

#include <iostream>
#include <vector>
using namespace std;

// Helper: copy a planar triangulation from MingCurve into DMWT.
static void applyPlanarTriangulation(MingCurve* myCurve, DMWT* myDMWT) {
	const Eigen::MatrixXi& triMat = myCurve->getPlanarTriangulation();
	int* triIndices = new int[triMat.rows() * 3];
	for (int i = 0; i < triMat.rows(); i++)
		for (int j = 0; j < 3; j++)
			triIndices[i*3 + j] = triMat(i, j);
	myDMWT->setTriangulation(triIndices, triMat.rows());
	myDMWT->buildList();
}

// interface 1
int delaunayRestrictedTriangulation(
	const double* inCurve, const int inNum,
	std::vector<double>& outPoints,
	std::vector<double>& outFaces,
	float* weights,
	bool dosmooth, int subd, int laps)
{
	const bool withNorm = false;
	const bool isdmwt   = true;
	const bool isdot1   = false;

	float weightTri   = weights[0];
	float weightEdge  = weights[1];
	float weightBiTri = weights[2];
	float weightTriBd = weights[3];
	float weightWorst = weights[4];
	const int limit   = 1000000;

	try {
		MingCurve* myCurve = new MingCurve(inCurve, inNum, limit, withNorm);
		if (!myCurve->edgeProtect(isdmwt)) {
			delete myCurve;
			Geex::cs_log() << "MWT: (0) bad input, not able to protect curve" << endl;
			return 0;
		}
		int     ptn      = myCurve->getNumOfPoints();
		double* pts      = myCurve->getPoints();
		double* deGenPts = myCurve->getDeGenPoints();

		DMWT* myDMWT = new DMWT(ptn, pts, deGenPts, myCurve->isDeGen);
		myDMWT->setWeights(weightTri, weightEdge, weightBiTri, weightTriBd, weightWorst);
		myDMWT->setDot(isdot1);

		if (myCurve->hasPlanarTriangulation())
			applyPlanarTriangulation(myCurve, myDMWT);
		else
			myDMWT->preprocess();

		if (!myDMWT->start()) {
			delete myDMWT;
			delete myCurve;
			Geex::cs_log() << "MWT: (0) no solution case" << endl;
			return 0;
		}
		std::vector<float> norms;  // discarded — interface 1 has no normals
		myDMWT->getResult(outFaces, outPoints, norms, dosmooth, subd, laps);
		delete myDMWT;
		delete myCurve;
	} catch (const std::exception& e) {
		Geex::cs_err() << "MWT: Exception: " << e.what() << endl;
		return 0;
	} catch (...) {
		Geex::cs_err() << "MWT: Unknown error" << endl;
		return 0;
	}
	return 1;
}

// interface 2
int delaunayRestrictedTriangulation(
	const double* inCurve, const float* inNorm, const int inNum,
	std::vector<double>& outPoints, std::vector<float>& outNorms,
	std::vector<double>& outFaces,
	float* weights,
	bool dosmooth, int subd, int laps)
{
	const bool withNorm = true;
	const bool isdmwt   = true;
	const bool isdot1   = false;

	float weightTri   = weights[0];
	float weightEdge  = weights[1];
	float weightBiTri = weights[2];
	float weightTriBd = weights[3];
	float weightWorst = weights[4];
	const int limit   = 1000000;

	try {
		MingCurve* myCurve = new MingCurve(inCurve, inNorm, inNum, limit, withNorm);
		if (!myCurve->edgeProtect(isdmwt)) {
			delete myCurve;
			Geex::cs_log() << "MWT: (0) bad input, not able to protect curve" << endl;
			return 0;
		}
		int     ptn      = myCurve->getNumOfPoints();
		double* pts      = myCurve->getPoints();
		double* deGenPts = myCurve->getDeGenPoints();
		float*  norms    = myCurve->getNormal();

		DMWT* myDMWT = new DMWT(ptn, pts, deGenPts, norms, myCurve->isDeGen);
		myDMWT->setWeights(weightTri, weightEdge, weightBiTri, weightTriBd, weightWorst);
		myDMWT->setDot(isdot1);

		if (myCurve->hasPlanarTriangulation())
			applyPlanarTriangulation(myCurve, myDMWT);
		else
			myDMWT->preprocess();

		if (!myDMWT->start()) {
			delete myDMWT;
			delete myCurve;
			Geex::cs_log() << "MWT: (0) no solution case" << endl;
			return 0;
		}
		myDMWT->getResult(outFaces, outPoints, outNorms, dosmooth, subd, laps);
		delete myDMWT;
		delete myCurve;
	} catch (const std::exception& e) {
		Geex::cs_err() << "MWT: Exception: " << e.what() << endl;
		return 0;
	} catch (...) {
		Geex::cs_err() << "MWT: Unknown error" << endl;
		return 0;
	}
	return 1;
}

// interface 3
int delaunayRestrictedTriangulation(
	const double* inCurve, const float* inNorm, const int inNum,
	std::vector<double>& outPoints, std::vector<float>& outNorms,
	std::vector<int>& outFaces,
	float* weights)
{
	const bool withNorm = true;
	const bool isdmwt   = true;
	const bool isdot1   = false;

	float weightTri   = weights[0];
	float weightEdge  = weights[1];
	float weightBiTri = weights[2];
	float weightTriBd = weights[3];
	float weightWorst = weights[4];
	const int limit   = 1000000;

	try {
		MingCurve* myCurve = new MingCurve(inCurve, inNorm, inNum, limit, withNorm);
		if (!myCurve->edgeProtect(isdmwt)) {
			delete myCurve;
			Geex::cs_log() << "MWT: (0) bad input, not able to protect curve" << endl;
			return 0;
		}
		int     ptn      = myCurve->getNumOfPoints();
		double* pts      = myCurve->getPoints();
		double* deGenPts = myCurve->getDeGenPoints();
		float*  norms    = myCurve->getNormal();

		DMWT* myDMWT = new DMWT(ptn, pts, deGenPts, norms, myCurve->isDeGen);
		myDMWT->setWeights(weightTri, weightEdge, weightBiTri, weightTriBd, weightWorst);
		myDMWT->setDot(isdot1);

		if (myCurve->hasPlanarTriangulation())
			applyPlanarTriangulation(myCurve, myDMWT);
		else
			myDMWT->preprocess();

		if (!myDMWT->start()) {
			delete myDMWT;
			delete myCurve;
			Geex::cs_log() << "MWT: (0) no solution case" << endl;
			return 0;
		}
		myDMWT->getResult(outFaces, outPoints, outNorms);
		delete myDMWT;
		delete myCurve;
	} catch (const std::exception& e) {
		Geex::cs_err() << "MWT: Exception: " << e.what() << endl;
		return 0;
	} catch (...) {
		Geex::cs_err() << "MWT: Unknown error" << endl;
		return 0;
	}
	return 1;
}