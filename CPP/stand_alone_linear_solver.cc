
#include <iostream>
#include <iomanip>

// Get header (and body) of actual matrix and linear solver
#include "stand_alone_linear_solver.h"


using namespace std;

//===================================================================
/// Matrix test code
//===================================================================
int main()
{

 // initialize random seed:
 srand (time(NULL));

 // Choose size of system
 unsigned n=100;
 // Matrix and exact solution
 DenseDoubleMatrix A(n,n,0.0);
 std::vector<double> x_exact(n,0.0);
 for (unsigned i=0;i<n;i++)
  {
  // Create exact solution
   x_exact[i]=double(i);
   for (unsigned j=0;j<n;j++)
    {
     A(i,j)=-50.0+100*double(rand()%100);
    }
  }
 std::cout << "\nHere's the matrix: \n";
 A.output(cout);
 
// Soln and rhs
 vector<double> x(n,0.0);
 vector<double> b(n,0.0);
 
// Cook up RHS for exact solution
 A.multiply(x_exact,b);

// Solve the bastard
 A.solve(b,x);

// Check result:
 std::cout << "\n\ni\tx[i]\t\t x_exact[i]\tError\n";
 std::cout << "--------------------------------------------------\n";
 double error=0.0;
 for (unsigned i=0;i<n;i++)
  {
   double err=fabs(x[i]-x_exact[i]);
   std::cout << i << "\t" << std::setprecision(5) 
             << std::setw(10) << x[i] << "\t\t" << x_exact[i] << "\t" << err 
             << std::endl;
   error+=err;
  }
 std::cout << "\nGlobal error: " << error << std::endl;
 
}
