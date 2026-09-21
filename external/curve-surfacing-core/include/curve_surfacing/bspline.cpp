#pragma once
#include "curve_surfacing_inline.h"

// ── implementation includes ──
#include "bspline.h"
#include <vector>



//use deboor algorithm to compute spline curves

#include <vector>


//methods to manage structs. ENCAPSULATION
CS_INLINE void BSplineCurve::InitializeBSplineCurve()
{
	Xcoeff=NULL;
	Ycoeff=NULL;
	Zcoeff=NULL;
	knot=NULL;
	Xpoints=NULL;
	Ypoints=NULL;
	Zpoints=NULL;
	weights=NULL;
}

void BSplineCurve::CleanBSplineCurve(){
	if(Xcoeff!=NULL)
		delete[] Xcoeff;
	if(Ycoeff!=NULL)
		delete[] Ycoeff;
	if(Zcoeff!=NULL)
		delete[] Zcoeff;
	if(weights!=NULL)
		delete[] weights;
	if(knot!=NULL)
		delete[] knot;
	if(Xpoints!=NULL)
		delete[] Xpoints;
	if(Ypoints!=NULL)
		delete[] Ypoints;
	if(Zpoints!=NULL)
		delete[] Zpoints;
	point_num=0;
	InitializeBSplineCurve();
}

void BSplineCurve::ComputeBSplineCurve(){
	bspl_to_points(degree,l,Xcoeff,Ycoeff,Zcoeff,weights,knot,dense,Xpoints,Ypoints,Zpoints,&point_num);
}

int BSplineCurve::DEFAULT_DENSE=20;
int BSplineCurve::MAX_CONTROL_POINTS=30;

CS_INLINE BSplineCurve::BSplineCurve() : weights(NULL)
{
	this->InitializeBSplineCurve();
}

BSplineCurve::~BSplineCurve(){
	this->CleanBSplineCurve();
}

void bspl_to_points(int degree,int l,double xcoeff[],double ycoeff[],double zcoeff[],double weights[],double knot[],int dense,double xpoints[],double ypoints[],double zpoints[],int *point_num)
{
	int i,ii;
	double u;
	double denominator;

	*point_num=0;

	for(i=degree;i<l+degree;i++)
	{
		if(knot[i+1]>knot[i])  /*skip zero length intervals*/
		{
			for(ii=0;ii<dense;ii++)
			{
				u=knot[i]+ii*(knot[i+1]-knot[i])/dense;
				
				denominator=deboor(degree,weights,knot,u,i);
				xpoints[*point_num]=deboor(degree,xcoeff,knot,u,i)/denominator;
				ypoints[*point_num]=deboor(degree,ycoeff,knot,u,i)/denominator;
				zpoints[*point_num]=deboor(degree,zcoeff,knot,u,i)/denominator;
				*point_num=(*point_num)+1;

			}
			if(i==l+degree-1)
			{
				u=knot[i+1];
				denominator=deboor(degree,weights,knot,u,i);
				xpoints[*point_num]=deboor(degree,xcoeff,knot,u,i)/denominator;
				ypoints[*point_num]=deboor(degree,ycoeff,knot,u,i)/denominator;
				zpoints[*point_num]=deboor(degree,zcoeff,knot,u,i)/denominator;
				*point_num=(*point_num)+1;

			}
		}
	}
}

double deboor(int degree,double coeff[],double knot[],double u,int i)
{
	int k,j;
	double t1,t2;
	std::vector<double> coeffa (i+1); //auxiliary array
	for(j=i-degree;j<=i;j++) coeffa[j]=coeff[j]; //save data

	for(k=1;k<=degree;k++)
		for(j=i;j>=i-degree+k;j--)
		{
			t1=(knot[j+degree-k+1]-u)/(knot[j+degree-k+1]-knot[j]);
			t2=1.0-t1;
			coeffa[j]=t1*coeffa[j-1]+t2*coeffa[j];
		}
		
	return coeffa[i];
}

void knot_insertion(int &n_points,int degree,int index,double value,double coeff[],double knot[])
{
	int start=index-degree;
	double newcoeff[30];

	int i;
	for(i=start+1;i<=index;i++)
	{
		double ai=value-knot[i];
		if(ai!=0)
			ai/=(knot[i+degree]-knot[i]);
		newcoeff[i-start-1]=(1-ai)*coeff[i-1]+ai*coeff[i];
	}

	for(i=n_points-1;i>=index;i--)
		coeff[i+1]=coeff[i];
	for(i=start+1;i<=index;i++)
		coeff[i]=newcoeff[i-start-1];

	n_points++;
}

CS_INLINE void BSplineCurve::KnotInsertion(int index, double value)
{
if(index<0)
{
	if(value<0) // for test
	index=degree;
	else{
	int i;
	for(i=degree;i<l+degree;i++)
		if(knot[i]<=value&&value<knot[i+1])
			break;
	index=i;
	}
}
if(value<0)
value=knot[index]+(knot[index+1]-knot[index])/2;


//copyed from knot_insertion() and edited
	int start=index-degree;
	double newcoeffX[30],newcoeffY[30],newcoeffZ[30],newweight[30];

	int i;
	for(i=start+1;i<=index;i++)
	{
		double ai=value-knot[i];
		if(ai!=0)
			ai/=(knot[i+degree]-knot[i]);
		newcoeffX[i-start-1]=(1-ai)*Xcoeff[i-1]+ai*Xcoeff[i];
		newcoeffY[i-start-1]=(1-ai)*Ycoeff[i-1]+ai*Ycoeff[i];
		newcoeffZ[i-start-1]=(1-ai)*Zcoeff[i-1]+ai*Zcoeff[i];
		newweight[i-start-1]=(1-ai)*weights[i-1]+ai*weights[i];
	}

	for(i=n_points-1;i>=index;i--)
	{Xcoeff[i+1]=Xcoeff[i];
	Ycoeff[i+1]=Ycoeff[i];
	Zcoeff[i+1]=Zcoeff[i];
	weights[i+1]=weights[i];
	}
	for(i=start+1;i<=index;i++)
	{	Xcoeff[i]=newcoeffX[i-start-1];
		Ycoeff[i]=newcoeffY[i-start-1];
		Zcoeff[i]=newcoeffZ[i-start-1];
		weights[i]=newweight[i-start-1];
	}

	for(i=n_points+degree;i>index;i--)
		knot[i+1]=knot[i];
	knot[i+1]=value;

	n_points++;
	l=n_points-degree;
}

CS_INLINE void BSplineCurve::MyCurveSplit(int index, double value)
{
	int i;
for(i=0;i<degree+1;i++)
	KnotInsertion(index+i,value);

n_points=n_points-index-1;
l=n_points-degree;
for(i=0;i<=n_points+degree;i++)
{
	knot[i]=knot[i+index+1];
}
for(i=0;i<n_points;i++)
{
	Xcoeff[i]=Xcoeff[i+index+1];
	Ycoeff[i]=Ycoeff[i+index+1];
	Zcoeff[i]=Zcoeff[i+index+1];
}

}



//compute the derivate of the curve, without the array of sampled points, since the derivatie is usually 
//not shown.

CS_INLINE BSplineCurve* BSplineCurve::DerivativeOrderOne()
{
	BSplineCurve *derivative=new BSplineCurve();
	derivative->degree=this->degree-1;
	derivative->dense=dense;
	derivative->n_points=n_points-1;
	derivative->l=derivative->n_points-derivative->degree;

	derivative->knot=new double[derivative->n_points+derivative->degree+1];
	for(int i=1;i<n_points+degree;i++)
		derivative->knot[i-1]=knot[i];

	

	derivative->Xcoeff=new double[derivative->n_points];
	derivative->Ycoeff=new double[derivative->n_points];
	derivative->Zcoeff=new double[derivative->n_points];
	derivative->weights=new double[derivative->n_points];
	for(int i=0;i<derivative->n_points;i++)
	{
		derivative->Xcoeff[i]=degree/(knot[i+degree+1]-knot[i+1])*(Xcoeff[i+1]-Xcoeff[i]);
		derivative->Ycoeff[i]=degree/(knot[i+degree+1]-knot[i+1])*(Ycoeff[i+1]-Ycoeff[i]);
		derivative->Zcoeff[i]=degree/(knot[i+degree+1]-knot[i+1])*(Zcoeff[i+1]-Zcoeff[i]);
		derivative->weights[i]=degree/(knot[i+degree+1]-knot[i+1])*(weights[i+1]-weights[i]);
	}
	
	return derivative;
}

double inline point_project_to_linesegment(double ptx, double pty, double ptz, double s0x, double s0y, double s0z, double s1x, double s1y, double s1z, double* dist){
	double s0s1_2 = (s0x-s1x)*(s0x-s1x) + (s0y-s1y)*(s0y-s1y) + (s0z-s1z)*(s0z-s1z);
	double ps0_s0s1 = (s0x-ptx)*(s1x-s0x) + (s0y-pty)*(s1y-s0y) + (s0z-ptz)*(s1z-s0z);
	double ps0_2 = (s0x-ptx)*(s0x-ptx) + (s0y-pty)*(s0y-pty) + (s0z-ptz)*(s0z-ptz);

	double t = -ps0_s0s1/s0s1_2;
	if (t < 0)
	{
		*dist = sqrt(ps0_2);
		return 0;
	}
	else if (t > 1)
	{
		double ps1_2 = (s1x-ptx)*(s1x-ptx) + (s1y-pty)*(s1y-pty) + (s1z-ptz)*(s1z-ptz);
		*dist = sqrt(ps1_2);
		return 1;
	}
	else{
		*dist = sqrt( abs(- t*t*s0s1_2 + ps0_2) ); //sqrt( t*t*s0s1_2 + 2*ps0_s0s1*t + ps0_2 );
		return t;
	}
}

// input the point to project, return the parameter of the projection
CS_INLINE double BSplineCurve::PointProjection(double x_point, double y_point, double z_point,double *x_proj,double *y_proj,double *z_proj)
{
	if (degree == 0)
	{
		return -1;
	}

	if(degree == 1){
		double min_dist = 1e10;
		double min_d_para = -1;
		double min_d_px, min_d_py, min_d_pz;
		for (int i = 0; i < n_points-1; i++)
		{
			double dist;
			double para = point_project_to_linesegment(x_point, y_point, z_point, Xcoeff[i], Ycoeff[i], Zcoeff[i], Xcoeff[i+1], Ycoeff[i+1], Zcoeff[i+1], &dist);

			if (dist < min_dist)
			{
				min_dist = dist;
				min_d_para = knot[i+1]*(1-para) + knot[i+2]*para;
				min_d_px = Xcoeff[i]*(1-para) + Xcoeff[i+1]*para;
				min_d_py = Ycoeff[i]*(1-para) + Ycoeff[i+1]*para;
				min_d_pz = Zcoeff[i]*(1-para) + Zcoeff[i+1]*para;
			}
		}

		*x_proj = min_d_px;
		*y_proj = min_d_py;
		*z_proj = min_d_pz;

		return min_d_para;
	}


	int point_count=0;
	int steps=10;
	double* auxXPoint=new double[steps*this->l];
	double* auxYPoint=new double[steps*this->l];
	double* auxZPoint=new double[steps*this->l];
	double* auxUArray=new double[steps*this->l];
	int i,ii;
	double u;
	for(i=this->degree;i<this->l+this->degree;i++)
	{
		if(this->knot[i+1]>this->knot[i])  /*skip zero length intervals*/
			for(ii=0;ii<steps;ii++)
			{
				u=knot[i]+ii*(knot[i+1]-knot[i])/steps;
				auxUArray[point_count]=u;
				double weight;
				weight=deboor(degree,weights,knot,u,i);
				auxXPoint[point_count]=deboor(degree,Xcoeff,knot,u,i)/weight;
				auxYPoint[point_count]=deboor(degree,Ycoeff,knot,u,i)/weight;
				auxZPoint[point_count]=deboor(degree,Zcoeff,knot,u,i)/weight;
				point_count=point_count+1;
			}
	}

	double distance_square=99999;
	int index=-1;
	for(i=0;i<point_count;i++)
	{
		double dist=(auxXPoint[i]-x_point)*(auxXPoint[i]-x_point)+(auxYPoint[i]-y_point)*(auxYPoint[i]-y_point)+
			(auxZPoint[i]-z_point)*(auxZPoint[i]-z_point);
		if(dist<distance_square)
		{distance_square=dist;
		index=i;
		}
	}

	//auxUArray[index] is the start value of newton iteration
	//to iterate, derivatives of first and second order are needed
	BSplineCurve *first_order=DerivativeOrderOne();
	BSplineCurve *second_order=first_order->DerivativeOrderOne();


	//thresholds of euclidean distance and cos zero are needed
	double thr_euc=0.01;
	double thr_cos=0.01;

	//ends of the knot vector are also needed
	double a=knot[degree];
	double b=knot[degree+l];

	//the iteration process
	do{
		//distance check
		double dist;
		dist=(auxXPoint[index]-x_point)*(auxXPoint[index]-x_point)+(auxYPoint[index]-y_point)*(auxYPoint[index]-y_point)+
			(auxZPoint[index]-z_point)*(auxZPoint[index]-z_point);
		dist=sqrt(dist);
		if(dist<=thr_euc)
			break;

		//cos zero check
		double first_devs[3];
		double second_devs[3];
		this->NURBSDerivativeIII(first_order,second_order,auxUArray[index],first_devs,second_devs);

		double numerator=first_devs[0]*(auxXPoint[index]-x_point)+first_devs[1]*(auxYPoint[index]-y_point)+
			first_devs[2]*(auxZPoint[index]-z_point);
		double abs_numerator=abs(numerator);
		double mode_deriv_first=first_devs[0]*first_devs[0]+first_devs[1]*first_devs[1]+first_devs[2]*first_devs[2];
		double denominator=dist*sqrt(mode_deriv_first);
		if(abs_numerator/denominator<=thr_cos)
			break;

		//since both the above conditions are not satiffied, iterate
		double aux_u=auxUArray[index];
		
		denominator=mode_deriv_first+second_devs[0]*(auxXPoint[index]-x_point)+second_devs[1]*(auxYPoint[index]-y_point)+
			second_devs[2]*(auxZPoint[index]-z_point);
		auxUArray[index]=auxUArray[index]-numerator/denominator;

		//bounded
		if(auxUArray[index]<a)
			auxUArray[index]=a;
		else if(auxUArray[index]>b)
			auxUArray[index]=b;

		double weight;
		GetValueAt(auxUArray[index],&auxXPoint[index],&auxYPoint[index],&auxZPoint[index],&weight);
		auxXPoint[index]/=weight;
		auxYPoint[index]/=weight;
		auxZPoint[index]/=weight;

		//check if para doesn't change significantly
		if(sqrt(mode_deriv_first)*abs(auxUArray[index]-aux_u)<=thr_euc)
			break;

		//continue iterate
		
	}while(1);

	*x_proj=auxXPoint[index];
	*y_proj=auxYPoint[index];
	*z_proj=auxZPoint[index];
	u=auxUArray[index];

	//free all the temporary spaces
	delete[] auxXPoint;
	delete[] auxYPoint;
	delete[] auxZPoint;
	delete[] auxUArray;
	delete first_order;
	delete second_order;

	return u;
}

// input the parameter value, output the coordinates of the respective point
CS_INLINE void BSplineCurve::GetValueAt(double doublePara, double* pXValue, double* pYValue, double* pZValue,double *pWeight)
{
	int i;
	if(doublePara==knot[l+degree])
		i=l+degree-1;
	else
	{for(i=degree;i<l+degree;i++)
		if(knot[i]<=doublePara&&doublePara<knot[i+1])
			break;
	}
	
	*pWeight=deboor(degree,weights,knot,doublePara,i);
	*pXValue=deboor(degree,Xcoeff,knot,doublePara,i);
	*pYValue=deboor(degree,Ycoeff,knot,doublePara,i);
	*pZValue=deboor(degree,Zcoeff,knot,doublePara,i);
	
}

// determines if a point is on the curve by first checking the convex property and then the distance from the point to the curve
CS_INLINE bool BSplineCurve::IsOnCurve(double x_value, double y_value, double z_value, double threshold)
{
	//convex property check
	if(!ConvexCheck(x_value,y_value,z_value))
		return false;

	//distance check
	double x_proj;
	double y_proj;
	double z_proj;
	double u=PointProjection(x_value,y_value,z_value,&x_proj,&y_proj,&z_proj);
	double dist=sqrt((x_value-x_proj)*(x_value-x_proj)+(y_value-y_proj)*(y_value-y_proj)
		+(z_value-z_proj)*(z_value-z_proj));
	if(dist<=threshold)
		return true;
	return false;
}

//since the routine is only used against two dimensional curves(curves on a plane), the convex hull is easy to determine
CS_INLINE bool BSplineCurve::ConvexCheck(double x_value, double y_value, double z_value)
{
	return true;
}

// split itself into two curves, return the first, and itself changes to the second part
CS_INLINE BSplineCurve* BSplineCurve::CurveSplitToT(double value)
{
	int i;
	for(i=degree;i<l+degree;i++)
		if(knot[i]<=value&&value<knot[i+1])
			break;

	int index=i;	
	for(i=0;i<degree+1;i++)
	KnotInsertion(index+i,value);
	
	BSplineCurve* firstpart=new BSplineCurve();
	firstpart->degree=degree;
	firstpart->n_points=index+1;
	firstpart->l=firstpart->n_points-firstpart->degree;
	firstpart->dense=dense;
	firstpart->knot=new double[firstpart->MAX_CONTROL_POINTS+firstpart->degree+1];
	firstpart->Xcoeff=new double[firstpart->n_points];
	firstpart->Ycoeff=new double[firstpart->n_points];
	firstpart->Zcoeff=new double[firstpart->n_points];
	firstpart->weights=new double[firstpart->n_points];
	for(i=0;i<firstpart->n_points;i++)
	{
		firstpart->Xcoeff[i]=Xcoeff[i];
		firstpart->Ycoeff[i]=Ycoeff[i];
		firstpart->Zcoeff[i]=Zcoeff[i];
		firstpart->weights[i]=weights[i];
	}
	for(i=0;i<firstpart->n_points+firstpart->degree+1;i++)
	{
		firstpart->knot[i]=knot[i];
	}
	firstpart->Xpoints=new double[(firstpart->MAX_CONTROL_POINTS-firstpart->degree)*firstpart->dense+1];
	firstpart->Ypoints=new double[(firstpart->MAX_CONTROL_POINTS-firstpart->degree)*firstpart->dense+1];
	firstpart->Zpoints=new double[(firstpart->MAX_CONTROL_POINTS-firstpart->degree)*firstpart->dense+1];
	firstpart->ComputeBSplineCurve();
	

	n_points=n_points-index-1;
	l=n_points-degree;
	for(i=0;i<=n_points+degree;i++)
	{
		knot[i]=knot[i+index+1];
	}
	for(i=0;i<n_points;i++)
	{
		Xcoeff[i]=Xcoeff[i+index+1];
		Ycoeff[i]=Ycoeff[i+index+1];
		Zcoeff[i]=Zcoeff[i+index+1];
		weights[i]=weights[i+index+1];
	}
	this->ComputeBSplineCurve();

	return firstpart;
}

CS_INLINE void BSplineCurve::NURBSDerivativeIII(BSplineCurve* first_dev,BSplineCurve* second_dev, double u, double* dev_i,double* dev_ii)
{
	double x,y,z,w;
	first_dev->GetValueAt(u,&x,&y,&z,&w);
	double cx,cy,cz,cw;
	GetValueAt(u,&cx,&cy,&cz,&cw);
	cx/=cw;
	cy/=cw;
	cz/=cw;

	dev_i[0]=(x-w*cx)/cw;
	dev_i[1]=(y-w*cy)/cw;
	dev_i[2]=(z-w*cz)/cw;

	double sx,sy,sz,sw;
	second_dev->GetValueAt(u,&sx,&sy,&sz,&sw);
	dev_ii[0]=(sx-(2*w*dev_i[0]+sw*cx))/cw;
	dev_ii[1]=(sy-(2*w*dev_i[1]+sw*cy))/cw;
	dev_ii[2]=(sz-(2*w*dev_i[2]+sw*cz))/cw;

}

// compute basis functions at u
CS_INLINE void BSplineCurve::BasisFunc(int index, double u, double* N)
{
	if(index<0)
	{
		int i;
		if(u==knot[l+degree])
			i=l+degree-1;
		else
		{for(i=degree;i<l+degree;i++)
			if(knot[i]<=u&&u<knot[i+1])
				break;
		}
		index=i;
	}

	double *left=new double[degree],*right=new double[degree];

	N[0]=1;
	for(int j=1;j<=degree;j++)
	{
		left[j-1]=u-knot[index+1-j];
		right[j-1]=knot[index+j]-u;

		double saved=0.0;
		for(int r=0;r<j;r++)
		{
			double temp=N[r]/(left[j-r-1]+right[r]);
			N[r]=saved+temp*right[r];
			saved=left[j-r-1]*temp;
		}
		N[j]=saved;

	}

	delete[] left;
	delete[] right;
}


CS_INLINE void BSplineCurve::PolylineTangent(double para, double* pXValue, double* pYValue, double* pZValue){
	if (degree==1)
	{
		for (int i = 0; i < n_points-1; i++)
		{
			if (knot[i+1]<=para && para <=knot[i+2])
			{
				double intvl = knot[i+2]-knot[i+1];

				*pXValue = (Xcoeff[i+1]-Xcoeff[i])/intvl;
				*pYValue = (Ycoeff[i+1]-Ycoeff[i])/intvl;
				*pZValue = (Zcoeff[i+1]-Zcoeff[i])/intvl;
				return;
			}
		}
	}
}