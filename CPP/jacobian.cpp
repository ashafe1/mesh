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

int main() {


	vector<double> X(2);
	vector<double> res(2);
	vector<double> res_p(X.size() * res.size());
	double eps = 0.1;

	for (int i=0; i<X.size(); i++) 
	{
		X[i] = 1;
		res[i] = X[i] + 3;
		//cout << res[i] << endl;
	}

	for (int i=0; i<X.size(); i++)
	{
		for (int j=0; j<X.size(); j++)
		{	
			if (j == i) res_p[X.size()*i + j] = res[i] + eps;
			else res_p[X.size()*i + j] = res[i];
			//cout << res_p[3*i +j] << endl;
		}
	}

	double** jacobian = get_jacobian(res, res_p, X.size(), eps);

	for (int i=0; i<X.size(); i++)
	{
		for (int j=0; j<X.size(); j++)
		{
			cout << jacobian[i][j] << endl;
		}
	}
	
	return 0;
}