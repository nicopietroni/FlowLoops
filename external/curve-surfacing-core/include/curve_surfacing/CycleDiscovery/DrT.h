#ifndef _DRT_H
#define _DRT_H

// For CycleDiscovery Project
// use Sum of Area Metric Only: weight = {1,0,0,0}, use interface 1
// use Sum of Edge Metric Only: weight = {0,1,0,0}, use interface 1
// use Sum of Dihedral Metric Only: weight = {0,0,1,0}, use interface 1
// use Dihedral & Board Normal Metric together: weight = {0,0,1,1} (two of 1's), use interface 2

// return type:
// 0: error, could be (a) bad input (b) no solution from MWT (c) unknown error
// 1: correct results
// 2: degenerated case, input is a plain, no result returned

#include "DMWT.h"
#include "MingCurve.h"
#include <stdio.h>
#include <vector>

// interface 1: curve only, no normals, with optional smoothing
// outFaces:  flattened triangle vertex positions, size = numTris * 9
// outPoints: curve points after edge-protection, size = numPoints * 3
int delaunayRestrictedTriangulation(
	const double* inCurve, const int inNum,
	std::vector<double>& outPoints,
	std::vector<double>& outFaces,
	float* weights,
	bool dosmooth, int subd, int laps);

// interface 2: curve + normals, with optional smoothing
// outFaces:  flattened triangle vertex positions, size = numTris * 9
// outPoints: curve points after edge-protection, size = numPoints * 3
// outNorms:  per-point normals, size = numPoints * 3; empty if no normals
int delaunayRestrictedTriangulation(
	const double* inCurve, const float* inNorm, const int inNum,
	std::vector<double>& outPoints, std::vector<float>& outNorms,
	std::vector<double>& outFaces,
	float* weights,
	bool dosmooth, int subd, int laps);

// interface 3: curve + normals, no smoothing, returns face indices
// outFaces:  flattened triangle vertex indices, size = numTris * 3
// outPoints: curve points after edge-protection, size = numPoints * 3
// outNorms:  per-point normals, size = numPoints * 3; empty if no normals
int delaunayRestrictedTriangulation(
	const double* inCurve, const float* inNorm, const int inNum,
	std::vector<double>& outPoints, std::vector<float>& outNorms,
	std::vector<int>& outFaces,
	float* weights);


#ifndef CURVE_SURFACING_STATIC
#  include "DrT.cpp"
#endif

#endif