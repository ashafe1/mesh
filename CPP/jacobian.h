//================================================
/// Header file for function to find the jacobian
//================================================

#ifndef JACOBIAN_HEADER
#define JACOBIAN_HEADER

#include <iostream>
using namespace std;
#include <math.h>
#include <vector>

//================================================================
/// Function to get the Jacboian matrix using a finite difference,
/// given the arguments: f(x), f(x+e), dimension, epsilon.
//================================================================

double** get_jacobian(vector<double> residuals, vector<double> residuals_perturb, unsigned dim, double eps)
{	
	/// Determine size of the Jacobian
	unsigned rows = residuals.size();
	unsigned cols = dim;

	/// Build the empty Jacobian
	double** jacobian = new double*[rows];
	for (unsigned i=0; i<rows; i++)
	{
		jacobian[i] = new double[cols];
	}

	/// Fill the Jacobian

	unsigned n=0; /// residuals_perturb index

	for (unsigned i=0; i<rows; i++)
	{
		for (unsigned j=0; j<cols; j++)
		{
			jacobian[i][j] = (residuals_perturb[n] - residuals[i]) / eps;
			n++;
		}
	}

	return jacobian;
}

#endif