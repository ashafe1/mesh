#include <math.h>
#include "mvector.h"
#include <vector>
using namespace std;


// // abstract base class
// struct MFunction
// {
// // pure virtual function
// virtual MVector operator()(const double& K,
// 						   const MVector& y, const MVector& z, const double& r) = 0;
// };

// // derived function class from above abstract base class
// class ForceLeft : public MFunction
// {
// public:
// 	MVector operator()(const double& K, const MVector& y, const MVector& z, const double& r)
// 	{
// 		MVector temp(2);
// 		temp[0] = K*(sqrt((y[1] - z[1])*(y[1] - z[1]) + (y[2] - z[2])*(y[2] - z[2])) - r)*((z[1]-y[1])/(sqrt((y[1] - z[1])*(y[1] - z[1]) + (y[2] - z[2])*(y[2] - z[2]))));
// 		temp[1] = x*y[0] - y[1];
// 		return temp;
// 	}
// };

int main() {

	int n = 3;  // Lattice size (n x n)

	double K = 3; // Spring constant (stiffness)

	int H = n-1; // Original height

	vector< MVector > P(n*n, MVector(2)); // Array of initial position vectors for each node

	// Fix the nodes on the floor at y = 0
	for(int i = 0; i < n; i++) {
		P[i][1] = i;
		P[i][2] = 0;
	}

	// Fix the nodes at y = H
	for(int i = n*(n-1), j = 0; i < n*n; i++, j++) {
		P[i][1] = j;
		P[i][2] = H;
	}

	// "Free" nodes
	int x = 0;
	int y = 1;
	for(int i = n; i < n*(n-1); i++) {
		if(x > n-1) {
			x = 0;
			y += 1;
		}
		P[i][1] = x;
		P[i][2] = y;

		x += 1;
	}

	vector< vector<double> > PFree(n*(n-2), vector<double>(2));  // Array containing all free nodes

	// Assign the free nodes
	for(int i = 0, j = n; i <= n*(n-2)-1 , j < n*(n-1); i++, j++) {
		PFree[i][1] = P[j][1];
		PFree[i][2] = P[j][2];
		} 
	
	cout << "All nodes" << '\n';
	// Output check that the initial position vector is correct
	for(int i = 0; i < n*n; i++) cout << "(" <<P[i][1]<<", "<<P[i][2]<<")"<<'\n';

	cout << "Free nodes" << '\n';
	// Output check that the initial free node vector is correct
	for(int i = 0; i < n*(n-2); i++) cout << "(" <<PFree[i][1]<<", "<<PFree[i][2]<<")"<<'\n';
	
	// Control parameter (new height)
	int h = 5;

	// New fixed node positions
	vector< vector<double> > p(n*n, vector<double>(2)); // Array of new position vectors for each node

		// Fix the nodes at y = h
	for(int i = n*(n-1), j = 0; i < n*n; i++, j++) {
		p[i][1] = j;
		p[i][2] = h;
	}

	// Keep floor nodes fixed
	for(int i = 0; i < n; i++) p[i] = P[i];

	// Output to check new position vectors
	cout << "New position vectors" << '\n';
	for(int i = 0; i < n*n; i++) cout << "(" <<p[i][1]<<", "<<p[i][2]<<")"<<'\n';	

	// Define sets of neighbours for each free node
	vector< vector< vector<double> > > neighbour(n*(n-2) + 1, vector<double>(2));
	for(int i = n; i < n*(n-1); i++) {
		if(i % n == 0) 
		{
			
		}
	}

	return 0;

}