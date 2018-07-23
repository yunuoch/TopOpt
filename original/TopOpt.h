#pragma once
#include <vector>
#include <math.h>
#include<iostream>
#include <GL/freeglut.h>
#include<Eigen\Dense>
#include<Eigen\Sparse>
//#include<Eigen\Cholesky>
#include <iterator>
//#include <algorithm>

using namespace std;
using namespace Eigen;

class TopOpt
{
public:
	TopOpt();
	TopOpt(int x, int y, double v, double p, double rm);
	virtual ~TopOpt();
	int nelx, nely;
	double volfrac, penal, rmin;
	
	MatrixXd U, K, F;
	MatrixXd Ue, Ke, X, Xold, dc;

	void ke();
	void compute();
	void FEA();
	void OC();
	void Filter();

	void draw();

	double volume;
};

