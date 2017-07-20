//==============================================================
// Header for matrix class and associated dense linear solver
// Hacked around from oomph-lib
//==============================================================

//Include guards to prevent multiple inclusion of the header
#ifndef STAND_ALONE_LINEAR_SOLVER_HEADER
#define STAND_ALONE_LINEAR_SOLVER_HEADER



#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h> 
#include <vector> 
#include <math.h> 


//Forward definition of the linear solver class
 class LinearSolver; 

//Forward definition of the DenseLU class
class DenseLU;




/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////


//=============================================================================
/// \short Abstract base class for matrices of doubles -- adds 
/// abstract interfaces for solving, LU decomposition and 
/// multiplication by vectors.
//=============================================================================
class DoubleMatrixBase
{
  protected:

 //Pointer to a linear solver
 DenseLU* Linear_solver_pt;

  public:

 /// (Empty) constructor. 
 DoubleMatrixBase() : Linear_solver_pt(0) {}

 /// Broken copy constructor
 DoubleMatrixBase(const DoubleMatrixBase& matrix) 
  {
   abort();
   //BrokenCopy::broken_copy("DoubleMatrixBase");
  } 
 
 /// Broken assignment operator
 void operator=(const DoubleMatrixBase&) 
  {
   abort();
   //BrokenCopy::broken_assign("DoubleMatrixBase");
  }

  /// Return the number of rows of the matrix
 virtual unsigned long nrow() const=0;
 
 /// Return the number of columns of the matrix
 virtual unsigned long ncol() const=0;

 /// virtual (empty) destructor
 virtual ~DoubleMatrixBase() { }

 /// \short Round brackets to give access as a(i,j) for read only
 /// (we're not providing a general interface for component-wise write
 /// access since not all matrix formats allow efficient direct access!)
 virtual double operator()(const unsigned long &i, 
                           const unsigned long &j) const=0;
 
 /// Return a pointer to the linear solver object
 DenseLU* &linear_solver_pt() {return Linear_solver_pt;}

 /// Return a pointer to the linear solver object (const version)
 DenseLU* const &linear_solver_pt() const {return Linear_solver_pt;}

 /// \short Complete LU solve (replaces matrix by its LU decomposition
 /// and overwrites RHS with solution). The default should not need
 /// to be over-written
 void solve(std::vector<double> &rhs);

 /// \short Complete LU solve (Nothing gets overwritten!). The default should
 /// not need to be overwritten
 void solve(const std::vector<double> &rhs, std::vector<double> &soln);

 /// \short Multiply the matrix by the vector x: soln=Ax.
 virtual void multiply(const std::vector<double> &x, 
                       std::vector<double> &soln)const=0;

};



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////



//======================================================================
/// \short Class for dense matrices, storing all the values of the 
/// matrix as a pointer to a pointer with assorted output functions 
/// inherited from Matrix<T>. The curious recursive template pattern is
/// used here to pass the specific class to the base class so that
/// round bracket access can be inlined.
//======================================================================
template <class T>
class DenseMatrix 
{

  protected:

 /// Internal representation of matrix as a pointer to data
 T* Matrixdata;

 /// Number of rows
 unsigned long N;
 
 /// Number of columns
 unsigned long M;
 
  public:

 /// Empty constructor, simply assign the lengths N and M to 0
 DenseMatrix() : Matrixdata(0), N(0), M(0) {}

 /// Copy constructor: Deep copy!
 DenseMatrix(const DenseMatrix& source_matrix) 
  {
   // Set row and column lengths
   N=source_matrix.nrow();
   M=source_matrix.ncol();
   // Assign space for the data
   Matrixdata = new T[N*M];
   // Copy the data across from the other matrix
   for(unsigned long i=0;i<N;i++)
    {
     for(unsigned long j=0;j<M;j++)
      {
       Matrixdata[M*i+j] = source_matrix(i,j);
      }
    }
  }

 /// Copy assignment 
 DenseMatrix& operator=(const DenseMatrix& source_matrix) 
  {
   // Don't create a new matrix if the assignment is the identity
   if (this != & source_matrix)
    {
     // Check row and column length
     unsigned long n=source_matrix.nrow();
     unsigned long m=source_matrix.ncol();
     if ( (N!=n) || (M!=m) )
      {
       resize(n,m);
      }
     // Copy entries across from the other matrix
     for (unsigned long i=0;i<N;i++)
      {
       for (unsigned long j=0;j<M;j++)
        {
         (*this)(i,j) = source_matrix(i,j);
        }
      }
    }
   // Return reference to object itself (i.e. de-reference this pointer)
   return *this;
  }
 
 /// \short The access function that will be called by the read-write
 /// round-bracket operator.
 inline T& entry(const unsigned long &i, const unsigned long &j) 
  {
#ifdef RANGE_CHECKING
   this->range_check(i,j);
#endif   
   return Matrixdata[M*i+j];
  }

 /// \short The access function the will be called by the read-only 
 /// (const version) round-bracket operator.
 inline T get_entry(const unsigned long &i, 
                    const unsigned long &j) const
  {
#ifdef RANGE_CHECKING
   this->range_check(i,j);
#endif  
   return Matrixdata[M*i+j];
  }

 /// \short Round brackets to give access as a(i,j) for read only
 /// (we're not providing a general interface for component-wise write
 /// access since not all matrix formats allow efficient direct access!)
 /// The function uses the  MATRIX_TYPE template parameter to call the
 /// get_entry() function which must be defined in all derived classes
 /// that are to be fully instantiated.
 inline T operator()(const unsigned long &i, 
                     const unsigned long &j) const
  {
   //return static_cast<MATRIX_TYPE const *>(this)->get_entry(i,j);
   return this->get_entry(i,j);
  }

 /// Constructor to build a square n by n matrix
 DenseMatrix(const unsigned long &n);
 
 /// Constructor to build a matrix with n rows and m columns
 DenseMatrix(const unsigned long &n, const unsigned long &m);
 
 /// \short Constructor to build a matrix with n rows and m columns,
 /// with initial value initial_val
 DenseMatrix(const unsigned long &n, const unsigned long &m,
             const T &initial_val);
 
 /// Destructor, clean up the matrix data
 virtual ~DenseMatrix() {delete[] Matrixdata; Matrixdata=0;}
 
 /// Return the number of rows of the matrix
 inline unsigned long nrow() const {return N;} 

 /// Return the number of columns of the matrix
 inline  unsigned long ncol() const {return M;}
 
 /// Resize to a square nxn matrix;
 /// any values already present will be transfered
 void resize(const unsigned long &n) {resize(n,n);}
 
 /// \short Resize to a non-square n x m matrix;
 /// any values already present will be transfered
 void resize(const unsigned long &n, const unsigned long &m);
 
 /// \short Resize to a non-square n x m matrix and initialize the 
 /// new values to initial_value.
 void resize(const unsigned long &n, const unsigned long &m, 
             const T& initial_value);
 
 /// \short Initialize all values in the matrix to val.
 void initialise(const T& val)
  {for(unsigned long i=0;i<(N*M);++i) {Matrixdata[i] = val;}}
 
 /// Output function to print a matrix row-by-row to the stream outfile
 void output(std::ostream &outfile) const;
 
 /// Output function to print a matrix row-by-row to a file. Specify filename.
 void output(std::string filename) const;

 /// \short Indexed output function to print a matrix to the 
 /// stream outfile as i,j,a(i,j)
 void indexed_output(std::ostream &outfile) const;

 /// \short Indexed output function to print a matrix to a
 /// file as i,j,a(i,j). Specify filename.
 void indexed_output(std::string filename) const;

 /// \short Output the "bottom right" entry regardless of it being
 /// zero or not (this allows automatic detection of matrix size in
 /// e.g. matlab, python).
 void output_bottom_right_zero_helper(std::ostream &outfile) const;

 /// \short Indexed output function to print a matrix to the 
 /// stream outfile as i,j,a(i,j) for a(i,j)!=0 only.
 void sparse_indexed_output_helper(std::ostream &outfile) const;

};


///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


//=================================================================
/// \short Class of matrices containing doubles, and stored as a 
/// DenseMatrix<double>, but with solving functionality inherited
/// from the abstract DoubleMatrix class. 
//=================================================================
class DenseDoubleMatrix : public DoubleMatrixBase,
                          public DenseMatrix<double>
{

  public:

 /// Constructor, set the default linear solver
 DenseDoubleMatrix();
 
 /// Constructor to build a square n by n matrix.
 DenseDoubleMatrix(const unsigned long &n);
 
 /// Constructor to build a matrix with n rows and m columns.
 DenseDoubleMatrix(const unsigned long &n, const unsigned long &m);

 /// \short Constructor to build a matrix with n rows and m columns,
 /// with initial value initial_val
 DenseDoubleMatrix(const unsigned long &n, const unsigned long &m,
                   const double &initial_val);

 /// Broken copy constructor
 DenseDoubleMatrix(const DenseDoubleMatrix& matrix) : 
  DoubleMatrixBase(), DenseMatrix<double>()
  {
   abort();
   //BrokenCopy::broken_copy("DenseDoubleMatrix");
  } 
 
 /// Broken assignment operator
 void operator=(const DenseDoubleMatrix&) 
  {
   abort();
   //BrokenCopy::broken_assign("DenseDoubleMatrix");
  }

 
 /// Return the number of rows of the matrix
 inline unsigned long nrow() const {return DenseMatrix<double>::nrow();}
 
 /// Return the number of columns of the matrix
 inline unsigned long ncol() const {return DenseMatrix<double>::ncol();}

 /// \short Overload the const version of the round-bracket access operator
 /// for read-only access.
 inline double operator()(const unsigned long &i, 
                          const unsigned long &j) 
  const {return DenseMatrix<double>::get_entry(i,j);}
   
 /// \short Overload the non-const version of the round-bracket access
 /// operator for read-write access
 inline double& operator()(const unsigned long &i, const unsigned long &j) 
  {return DenseMatrix<double>::entry(i,j);}
 
 /// Destructor
 virtual ~DenseDoubleMatrix(); 
 
 /// \short LU decomposition using DenseLU (default linea solver)
 virtual void ludecompose();

 /// \short LU backsubstitution
 virtual void lubksub(std::vector<double> &rhs);

 /// \short Multiply the matrix by the vector x: soln=Ax.
 virtual void multiply(const std::vector<double> &x, 
                       std::vector<double> &soln) const;
};


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////



//=============================================================================
/// \short Dense LU decomposition-based solve of full assembled linear system.
/// VERY inefficient but useful to illustrate the principle.
/// Only suitable for use with Serial matrices and vectors.
/// This solver will only work with non-distributed matrices and vectors
/// (note: DenseDoubleMatrix is not distributable)
//============================================================================
 class DenseLU 
{
 ///The DenseDoubleMatrix class is a friend
 friend class DenseDoubleMatrix;
 
  public:

 /// Constructor, initialise storage
 DenseLU() : 
  Sign_of_determinant_of_matrix(0),
  Index(0), LU_factors(0)
  {}

 /// Broken copy constructor
 DenseLU(const DenseLU& dummy)
  { 
   abort();
   //BrokenCopy::broken_copy("DenseLU");
  } 
 
 /// Broken assignment operator
 void operator=(const DenseLU&) 
  {
   abort();
   //BrokenCopy::broken_assign("DenseLU");
  }

 /// Destructor, clean up the stored LU factors
 ~DenseLU() {clean_up_memory();}

 /// \short Linear-algebra-type solver: Takes pointer to a matrix and rhs 
 /// vector and returns the solution of the linear system.
 void solve(DoubleMatrixBase* const &matrix_pt,const std::vector<double> &rhs,
            std::vector<double> &result);

  protected:

 /// Perform the LU decomposition of the matrix
 void factorise(DoubleMatrixBase* const &matrix_pt);

 /// Do the backsubstitution step to solve the system LU result = rhs
 void backsub(const std::vector<double> &rhs,
              std::vector<double> &result);

 /// Clean up the stored LU factors
 void clean_up_memory();

 /// \short Sign of the determinant of the matrix 
 /// (obtained during the LU decomposition)
 int Sign_of_determinant_of_matrix;
 
  private:

 /// Pointer to storage for the index of permutations in the LU solve
 long* Index;
 
 /// Pointer to storage for the LU decomposition
 double* LU_factors;
};


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//============================================================================
/// Complete LU solve (Nothing gets overwritten!). This generic
/// version should never need to be overwritten
//============================================================================
void DoubleMatrixBase::solve(const std::vector<double> &rhs, 
                             std::vector<double> &soln)
{
 //Use the linear algebra interface to the linear solver
 Linear_solver_pt->solve(this,rhs,soln);
}

//============================================================================
/// Complete LU solve (overwrites RHS with solution). This is the
/// generic version which should not need to be over-written.
//============================================================================
void DoubleMatrixBase::solve(std::vector<double> &rhs)
{
 // Copy rhs vector into local storage so it doesn't get overwritten
 // if the linear solver decides to initialise the solution vector, say,
 // which it's quite entitled to do!
 std::vector<double> actual_rhs(rhs);

 //Use the linear algebra interface to the linear solver
 Linear_solver_pt->solve(this,actual_rhs,rhs);
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

//===============================================================
/// Constructor, set the default linear solver to be the DenseLU 
/// solver
//===============================================================
DenseDoubleMatrix::DenseDoubleMatrix(): DenseMatrix<double>()
{
 Linear_solver_pt = new DenseLU;
}

//==============================================================
/// Constructor to build a square n by n matrix.
/// Set the default linear solver to be DenseLU
//==============================================================
DenseDoubleMatrix::DenseDoubleMatrix(const unsigned long &n) : 
 DenseMatrix<double>(n)
{
 Linear_solver_pt = new DenseLU;
}
 

//=================================================================
/// Constructor to build a matrix with n rows and m columns.
/// Set the default linear solver to be DenseLU
//=================================================================
 DenseDoubleMatrix::DenseDoubleMatrix(const unsigned long &n, 
                                      const unsigned long &m) :
  DenseMatrix<double>(n,m)
{
 Linear_solver_pt = new DenseLU;
}

//=====================================================================
/// Constructor to build a matrix with n rows and m columns,
/// with initial value initial_val
/// Set the default linear solver to be DenseLU
//=====================================================================
DenseDoubleMatrix::DenseDoubleMatrix(const unsigned long &n, 
                                     const unsigned long &m,
                                     const double &initial_val) :
 DenseMatrix<double>(n,m,initial_val) 
{
 Linear_solver_pt = new DenseLU;
}

//=======================================================================
/// Destructor
//======================================================================
DenseDoubleMatrix::~DenseDoubleMatrix()
{
}

//============================================================================
/// LU decompose a matrix, by using the default linear solver
/// (DenseLU)
//============================================================================
void DenseDoubleMatrix::ludecompose()
{
 //Use the default (DenseLU) solver to ludecompose the matrix
 Linear_solver_pt->factorise(this);
}


//============================================================================
///  Back substitute an LU decomposed matrix.
//============================================================================
void DenseDoubleMatrix::lubksub(std::vector<double> &rhs)
{
 //Use the default (DenseLU) solver to perform the backsubstitution
 Linear_solver_pt->backsub(rhs,rhs);
}



//============================================================================
///  Multiply the matrix by the vector x: soln=Ax
//============================================================================
void DenseDoubleMatrix::multiply(const std::vector<double> &x, 
                                 std::vector<double> &soln) const
{

 // Multiply the matrix A, by the vector x 
 const double* x_pt = &x[0];
 double* soln_pt = &soln[0];
 for (unsigned long i=0;i<N;i++)
  {
   soln[i]=0.0;
   for (unsigned long j=0;j<M;j++)
    {
     soln_pt[i] += Matrixdata[M*i+j]*x_pt[j];
    }
  }
}




/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


//============================================================================
/// Constructor to build a square n by n matrix
//============================================================================
template<class T>
DenseMatrix<T>::DenseMatrix(const unsigned long &n)
{
 // Set row and column lengths
 N=n; M=n;
 // Assign space for the n rows
 Matrixdata = new T[n*n];
 //Initialise to zero if required
#ifdef OOMPH_INITIALISE_DENSE_MATRICES
 initialise(T(0));
#endif
}


//============================================================================
/// Constructor to build a matrix with n rows and m columns
//============================================================================
template<class T>
DenseMatrix<T>::DenseMatrix(const unsigned long &n, 
                            const unsigned long &m)
{
 // Set row and column lengths
 N=n; M=m;
 // Assign space for the n rows
 Matrixdata = new T[n*m];
#ifdef OOMPH_INITIALISE_DENSE_MATRICES
 initialise(T(0));
#endif
}

//============================================================================
/// \short Constructor to build a matrix with n rows and m columns,
/// with initial value initial_val
//============================================================================
template<class T>
DenseMatrix<T>::DenseMatrix(const unsigned long &n, const unsigned long &m,
                            const T &initial_val)
{
 // Set row and column lengths
 N=n; M=m;
 // Assign space for the n rows
 Matrixdata = new T[n*m];
 initialise(initial_val);
}


//============================================================================
/// \short Resize to a non-square n_row x m_col matrix,
/// where any values already present will be transfered.
//============================================================================
template<class T>
void DenseMatrix<T>::resize(const unsigned long &n,
                            const unsigned long &m)
{
 //If the sizes are the same, do nothing
 if((n==N) && (m==M)) {return;}
 // Store old sizes
 unsigned long n_old=N, m_old=M;
 // Reassign the sizes
 N=n; M=m;
 // Store double pointer to old matrix data
 T* temp_matrix = Matrixdata;

 // Re-create Matrixdata in new size
 Matrixdata = new T[n*m];
 //Initialise to zero
#ifdef OOMPH_INITIALISE_DENSE_MATRICES 
 initialise(T(0));
#endif

 // Transfer previously existing values
 unsigned long n_copy, m_copy;
 n_copy = std::min(n_old,n); m_copy = std::min(m_old,m);

 // If matrix has values, transfer them to new matrix
 // Loop over rows
 for(unsigned long i=0;i<n_copy;i++)
  {
   // Loop over columns
   for (unsigned long j=0;j<m_copy;j++)
    {
     // Transfer values from temp_matrix
     Matrixdata[m*i+j] = temp_matrix[m_old*i+j];
    }
  }
 
 // Now kill storage for old matrix
 delete[] temp_matrix;
}



//============================================================================
/// \short Resize to a non-square n_row x m_col matrix and initialize the 
/// new entries to specified value.
//============================================================================
template<class T>
void DenseMatrix<T>::resize(const unsigned long &n,
                            const unsigned long &m, 
                            const T &initial_value)
{
 //If the size is not changed, just return
 if((n==N) && (m==M)) {return;}
 // Store old sizes
 unsigned long n_old=N, m_old=M;
 // Reassign the sizes
 N=n; M=m;
 // Store double pointer to old matrix data
 T* temp_matrix = Matrixdata;
 // Re-create Matrixdata in new size
 Matrixdata = new T[n*m];
 // Assign initial value (will use the newly allocated data)
 initialise(initial_value);
 
 // Transfering values
 unsigned long n_copy, m_copy;
 n_copy = std::min(n_old,n); m_copy = std::min(m_old,m);
 // If matrix has values, transfer them to temp_matrix
 // Loop over rows
 for (unsigned long i=0;i<n_copy;i++)
  {
   // Loop over columns
   for (unsigned long j=0;j<m_copy;j++)
    {
     // Transfer values to temp_matrix
     Matrixdata[m*i+j] = temp_matrix[m_old*i+j];
    }
  }

 // Now kill storage for old matrix
 delete[] temp_matrix;
}



//============================================================================
/// Output function to print a matrix row-by-row to the stream outfile
//============================================================================
template<class T>
void DenseMatrix<T>::output(std::ostream &outfile) const
{
 //Loop over the rows
 for(unsigned i=0;i<N;i++)
  {
   //Loop over the columne
   for(unsigned j=0;j<M;j++)
    {
     outfile << (*this)(i,j) << " ";
    }
   //Put in a newline
   outfile << std::endl;
  }
}



//============================================================================
/// Output function to print a matrix row-by-row to a file. Specify filename.
//============================================================================
template<class T>
void DenseMatrix<T>::output(std::string filename) const
{
 // Open file
 std::ofstream some_file;
 some_file.open(filename.c_str());
   
 output(some_file);
 some_file.close();
}



//============================================================================
/// Indexed output as i,j,a(i,j)
//============================================================================
template<class T>
void DenseMatrix<T>::indexed_output(std::ostream &outfile) const
{
 //Loop over the rows
 for(unsigned i=0;i<N;i++)
  {
   //Loop over the columns
   for(unsigned j=0;j<M;j++)
    {
     outfile << i << " " << j << " " << (*this)(i,j) << std::endl;
    }
  }
}



//============================================================================
/// \short Indexed output function to print a matrix to a
/// file as i,j,a(i,j). Specify filename.
//============================================================================
template<class T>
void DenseMatrix<T>::indexed_output(std::string filename) const
{
 // Open file
 std::ofstream some_file;
 some_file.open(filename.c_str());
 indexed_output(some_file);
 some_file.close();
}


//============================================================================
///Output the "bottom right" entry regardless of it being
/// zero or not (this allows automatic detection of matrix size in
/// e.g. matlab, python). 
//============================================================================
template<class T>
void DenseMatrix<T>::output_bottom_right_zero_helper(std::ostream &outfile) const
{
  int last_row = this->N-1;
  int last_col = this->M-1;

  // Use this strange thingy because of the CRTP discussed above.
  T last_value = this->operator()(last_row, last_col);
  
  if(last_value == T(0))
   {
    outfile << last_row << " " << last_col << " " << T(0)
            << std::endl;
   }
}

//============================================================================
/// Sparse indexed output as i,j,a(i,j) for a(i,j)!=0 only.
//============================================================================
template<class T>
void DenseMatrix<T>::sparse_indexed_output_helper(std::ostream &outfile) const
{
 //Loop over the rows
 for(unsigned i=0;i<N;i++)
  {
   //Loop over the column
   for(unsigned j=0;j<M;j++)
    {
     if ((*this)(i,j)!=T(0))
      {
       outfile << i << " " << j << " " << (*this)(i,j) << std::endl;
      }
    }
  }
}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////



//=============================================================================
/// LU decompose the matrix.
/// WARNING: this class does not perform any PARANOID checks on the std::vectors - 
/// these are all performed in the solve(...) method.
//=============================================================================
void DenseLU::factorise(DoubleMatrixBase* const &matrix_pt)
{
 //Set the number of unknowns
 const unsigned long n = matrix_pt->nrow();
 
 //Small constant
 const double small_number=1.0e-20;

 //Vector scaling stores the implicit scaling of each row
 std::vector<double> scaling(n);

 //Integer to store the sign that must multiply the determinant as
 //a consequence of the row/column interchanges
 int signature = 1;

 //Loop over rows to get implicit scaling information
 for(unsigned long i=0;i<n;i++)
  {
   double largest_entry=0.0;
   for(unsigned long j=0;j<n;j++)
    {
     double tmp = fabs((*matrix_pt)(i,j));
     if(tmp > largest_entry) largest_entry = tmp;
    }
   if(largest_entry==0.0) 
    {
     abort();
     // throw OomphLibError("Singular Matrix",
     //                     OOMPH_CURRENT_FUNCTION,
     //                     OOMPH_EXCEPTION_LOCATION);
    }
   //Save the scaling
   scaling[i] = 1.0/largest_entry;
  }

 //Firsly, we shall delete any previous LU storage.
 //If the user calls this function twice without changing the matrix
 //then it is their own inefficiency, not ours (this time).
 clean_up_memory();

 //Allocate storage for the LU factors, the index and store
 //the number of unknowns
 LU_factors = new double[n*n];
 Index = new long[n];

 //Now we know that memory has been allocated, copy over
 //the matrix values
 unsigned count=0;
 for(unsigned long i=0;i<n;i++)
  {
   for(unsigned long j=0;j<n;j++)
    {
     LU_factors[count] = (*matrix_pt)(i,j);
     ++count;
    }
  }

 //Loop over columns
 for(unsigned long j=0;j<n;j++)
  {
   //Initialise imax
   unsigned long imax=0;

   for(unsigned long i=0;i<j;i++)
    {
     double sum = LU_factors[n*i+j];
     for(unsigned long k=0;k<i;k++) 
      {
       sum -= LU_factors[n*i+k]*LU_factors[n*k+j];
      }
     LU_factors[n*i+j] = sum;
    }

   //Initialise search for largest pivot element
   double largest_entry=0.0;
   for(unsigned long i=j;i<n;i++)
    {
     double sum = LU_factors[n*i+j];
     for(unsigned long k=0;k<j;k++) 
      {
       sum -= LU_factors[n*i+k]*LU_factors[n*k+j];
      }
     LU_factors[n*i+j] = sum;
     //Set temporary
     double tmp = scaling[i]*fabs(sum);
     if(tmp >= largest_entry)
      {
       largest_entry = tmp;
       imax = i;
      }
    }

   //Test to see if we need to interchange rows
   if(j != imax)
    {
     for(unsigned long k=0;k<n;k++)
      {
       double tmp = LU_factors[n*imax+k];
       LU_factors[n*imax+k] = LU_factors[n*j+k];
       LU_factors[n*j+k] = tmp;
      }
     //Change the parity of signature
     signature = -signature;

     //Interchange scale factor
     scaling[imax] = scaling[j];
    }
   
   //Set the index
   Index[j] = imax;
   if(LU_factors[n*j+j] == 0.0) 
    {
     LU_factors[n*j+j] = small_number;
    }
   //Divide by pivot element
   if(j != n-1)
    {
     double tmp = 1.0/LU_factors[n*j+j];
     for(unsigned long i=j+1;i<n;i++) 
      {
       LU_factors[n*i+j] *= tmp;
      }
    }
  
  } //End of loop over columns

 
 //Now multiply all the diagonal terms together to get the determinant
 //Note that we need to use the mantissa, exponent formulation to
 //avoid underflow errors
 double determinant_mantissa=1.0;
 int determinant_exponent = 0, iexp;
 for(unsigned i=0; i<n; i++)
  {
   //Multiply by the next diagonal entry's mantissa
   //and return the exponent
   determinant_mantissa *= frexp(LU_factors[n*i+i], &iexp);

   //Add the new exponent to the current exponent
   determinant_exponent += iexp;

   // normalise
   determinant_mantissa = frexp(determinant_mantissa,&iexp);
   determinant_exponent += iexp;
  }

 //If paranoid issue a warning that the matrix is near singular
// #ifdef PARANOID
//  int tiny_exponent = -60;
//  if(determinant_exponent < tiny_exponent)
//   {
//    std::ostringstream warning_stream;
//    warning_stream << "The determinant of the matrix is very close to zero.\n"
//                   << "It is " << determinant_mantissa << " x 2^" 
//                   << determinant_exponent << "\n";
//    warning_stream << "The results will depend on the exact details of the\n"
//                   << "floating point implementation ... just to let you know\n";
//    OomphLibWarning(warning_stream.str(),
//                    "DenseLU::factorise()",
//                    OOMPH_EXCEPTION_LOCATION);
//   }
// #endif

 //Integer to store the sign of the determinant
 int sign = 0;

 //Find the sign of the determinant
 if(determinant_mantissa > 0.0) {sign = 1;}
 if(determinant_mantissa < 0.0) {sign = -1;}
 
 //Multiply the sign by the signature
 sign *= signature;
 
 //Return the sign of the determinant
 Sign_of_determinant_of_matrix = sign;
 }


//=============================================================================
/// Do the backsubstitution for the DenseLU solver.
/// WARNING: this class does not perform any PARANOID checks on the vectors -
/// these are all performed in the solve(...) method. So, if you call backsub
/// directly, you have been warned...
//=============================================================================
void DenseLU::backsub(const std::vector<double> &rhs,
                      std::vector<double> &result)
{
 //Copy the rhs vector into the result vector
 const unsigned long n = rhs.size();
 for(unsigned long i=0;i<n;++i)
  {
   result[i] = rhs[i];
  }
 
 // Loop over all rows for forward substition
 unsigned long k=0;
 for(unsigned long i=0;i<n;i++)
  {
   unsigned long ip = Index[i];
   double sum = result[ip];
   result[ip] = result[i];
   if(k != 0)
    {
     for(unsigned long j=k-1;j<i;j++)
      {
       sum -= LU_factors[n*i+j]*result[j];
      }
    }
   else if(sum != 0.0)
    {
     k = i+1;
    }
   result[i] = sum;
  }
 
  //Now do the back substitution
  for (long i=long(n)-1;i>=0;i--)
   {
    double sum = result[i];
    for(long j=i+1;j<long(n);j++)
     {
      sum -= LU_factors[n*i+j]*result[j];
     }
    result[i] = sum/LU_factors[n*i+i];
   }
}


//=============================================================================
 /// \short Linear-algebra-type solver: Takes pointer to a matrix and rhs 
 /// vector and returns the solution of the linear system. 
//============================================================================
 void DenseLU::solve(DoubleMatrixBase* const &matrix_pt,
                     const std::vector<double> &rhs,
                     std::vector<double> &result)
{
 
 // factorise
 factorise(matrix_pt);
 
  // backsubstitute
 backsub(rhs,result);
}
 

//=============================================================================
/// Delete the storage that has been allocated for the LU factors, if
/// the matrix data is not itself being overwritten.
//=============================================================================
void DenseLU::clean_up_memory()
{
 //Clean up the LU factor storage, if it has been allocated
 //N.B. we don't need to check the index storage as well.
 if(LU_factors!=0)
  {
   //Delete the pointer to the LU factors
   delete[] LU_factors;
   //Null out the vector
   LU_factors = 0;
   //Delete the pointer to the Index
   delete[] Index;
   //Null out
   Index=0;
  }
}


#endif
