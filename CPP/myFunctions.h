#ifndef myFunctions
#define myFunctions
#include <vector>
using namespace std;

///======================
/// Function prototypes
///======================

/// Find the number of nodes given the spatial vector
unsigned num_nodes(const std::vector<unsigned> nX);

/// Find the number of springs given the spatial vector
unsigned num_springs(std::vector<unsigned> nX);

/// Find the number of springs given the spatial vector
// unsigned num_springs(unsigned nx, unsigned ny, unsigned nz);


/// Find the norm of the difference between two vectors (spring length)
double diff_norm(const std::vector<double> u, const std::vector<double> v);






#endif