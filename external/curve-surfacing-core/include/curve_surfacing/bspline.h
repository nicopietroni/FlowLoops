#ifndef _bspline
#define _bspline

#include <stdlib.h>
#include <math.h>

class BSplineCurve {
public:
	//input
	int degree;
	int n_points;
	int l;
	double* Xcoeff;
	double* Ycoeff;
	double* Zcoeff;
	double* knot;
	int dense;
	//output
	double* Xpoints;
	double* Ypoints;
	double* Zpoints;
	int point_num;

public:
	// a test method for splitting curve
	void MyCurveSplit(int index,double value);

	// a const num, the maximum number of control points the datastructure can hold
	static int MAX_CONTROL_POINTS;

	// insert a knot into the knot vector without changing the shape of the curve.
	//input: index is the start knot of the range in which u lies
	// value is the u to be inserted as the new knot
	void KnotInsertion(int index,double value);

	BSplineCurve();
	~BSplineCurve();

	void InitializeBSplineCurve();

	void CleanBSplineCurve();

	// if the characteristics of the curve is changed, call this method to recalculate
	// the sample points.
	void ComputeBSplineCurve();
public:
	// this is function for non-rational bsplines. it computes the first order derivative of 
	// the curve, which is a new curve of lower degree
	BSplineCurve* DerivativeOrderOne();
public:
	// input the point to project, return the parameter of the projection
	double PointProjection(double x_point, double y_point, double z_point,double *x_proj,double *y_proj,double *z_proj);
public:
	// input the parameter value, output the coordinates of the respective point
	void GetValueAt(double doublePara, double* pXValue, double* pYValue, double* pZValue,double *pWeight);
public:
	// determines if a point is on the curve by first checking the convex property and then the distance from the point to the curve
	bool IsOnCurve(double x_value, double y_value, double z_value, double threshold);
public:
	// checks if the point is in the covex polygon specified by all control points of the spline curve
	bool ConvexCheck(double x_value, double y_value, double z_value);
public:
	// split itself into two curves, return the first, and itself changes to the second part
	BSplineCurve* CurveSplitToT(double value);
public:
	//weights of each control point for a nurbs curve.
	double* weights;
public:
	//it computes the first and second order derivative of the current curve at the specified u value
	//input: first_dev is the first derivative of the nurbs curve without being projected, that is, in the glasman space
	//		second_dev is the second order derivative of ......
	//		u is the parameter value
	//output:	dev_i is the vector of first order derivative
	//			dev_ii is the vector of second order derivative
	void NURBSDerivativeIII(BSplineCurve* first_dev,BSplineCurve* second_dev, double u, double* dev_i,double* dev_ii);
public:
	//default dense, number of sample points of a knot interval
	static int DEFAULT_DENSE;
public:
	// compute basis functions at u
	//output is in array N
	void BasisFunc(int index, double u, double* N);

	//compute tangent at para when the curve is degree one polyline.
	void PolylineTangent(double para, double* pXValue, double* pYValue, double* pZValue);
};




/*generate points on B-spline curve.(one coordinate)
Input: degree: polynomial degree of each piece of curve
		l:     number of active intervals
		coeff: B-spline control points

		knot:  knot sequence:knot[0]...knot[l+2*degree]
		dense: how many points per segment
Output: points:output array with function values
		*point_num: how many points are generated.
*/
void bspl_to_points(int degree,int l,double xcoeff[],double ycoeff[],double zcoeff[],double weights[],double knot[],int dense,double xpoints[],double ypoints[],double zpoints[],int *point_num);
/*uses deboor algorithm to compute one coordinate on B-spline curve for param. value u in interval i~i+1;
Input:	degree:	polynomial degree of each piece of curve
		coeff:	B-spline control points
		knot:	knot sequence
		u:	evaluation abscissa
		i:	u's interval:u[i]<=u<u[i+1]
Output:	coordinate value.
*/
double deboor(int degree,double coeff[],double knot[],double u,int i);

/*uses boehm algorithm to do one knot insertion at a time.
Input:	n_points: num of control points
		degree: polynomial degree of each piece of curve
		index:	index of the knot vector where new knot is inserted
		value:	the value of the knot to be inserted
		coeff:	B-spline control points
		knot:	knot sequence:knot[0]...knot[l+2*degree]
Output:	updated control points(coeff) and knots.
REQUIREMENT: the vector of coeff and knot should leave enough place for the 
renewed version.
*/
void knot_insertion(int &n_points,int degree,int index,double value,double coeff[],double knot[]);



#ifndef CURVE_SURFACING_STATIC
#  include "bspline.cpp"
#endif

#endif