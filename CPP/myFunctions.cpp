#include "myFunctions.h"
#include <vector>
#include <math.h>
using namespace std;

///=======================
/// Function definitions
///=======================

unsigned num_nodes(const std::vector<unsigned> nX)
{
  unsigned nNodes = 1;
  unsigned dim = nX.size();

  for (int i=0; i<dim; i++)
  {
    nNodes *= nX[i];
  }
  return nNodes;
}

// Find the number of springs given the spatial vector
unsigned num_springs(vector<unsigned> nX)
{
	unsigned dim = nX.size();
	unsigned nx = nX[0];
	unsigned ny = nX[1];
	if (dim == 2)
	{
		return ny*(nx-1) + nx*(ny-1);
	}
	else if (dim == 3)
	{
		unsigned nz = nX[2];
		return (nx*ny)*(nz-1) + nz*(ny*(nx-1) + nx*(ny-1));
	}
}

unsigned num_springs(unsigned nx, unsigned ny, unsigned nz)
{
	return (nx*ny)*(nz-1) + nz*(ny*(nx-1) + nx*(ny-1));
}


double diff_norm(const std::vector<double> u, const std::vector<double> v)
{
	double temp = 0;
	unsigned dim = u.size();
	for(int i = 0; i< dim; i++) temp += ((u[i] - v[i])*(u[i] - v[i]));
	return sqrt(temp);
}

