#ifndef MVECTOR_H // the 'include guard'
#define MVECTOR_H // see C++ Primer Sec. 2.9.2
#include <iostream>
#include <vector>
#include <stdlib.h>



// Class that represents a mathematical vector
class MVector
{
public:
	// constructors
	MVector() {}
	explicit MVector(int n) : v(n) {}
	MVector(int n, double x) : v(n, x) {}

	// access element (lvalue)
	double &operator[](int index) { return v[index]; }

	// access element (rvalue)
	double operator[](int index) const { return v[index]; }

	int size() const { return v.size(); } // number of elements in vector

private:
	std::vector<double> v;
};



#endif


// Operator overload for "scalar * vector"
inline MVector operator*(const double& lhs, const MVector& rhs)
{
	MVector temp(rhs);
	for (int i=0; i<temp.size(); i++) temp[i]*=lhs;
	return temp;
}

// Operator overload for "vector * scalar"
inline MVector operator*(const MVector& lhs, const double& rhs)
{
	MVector temp(lhs);
	for (int i=0; i<temp.size(); i++) temp[i]*=rhs;
	return temp;
}

// Operator overload for "vector / scalar"
inline MVector operator/(const MVector& lhs, const double& rhs)
{
	MVector temp(lhs);
	for (int i=0; i<temp.size(); i++) temp[i]/=rhs;
	return temp;
}

// Operator overload for "vector + vector"
inline MVector operator+(const MVector& lhs, const MVector& rhs)
{
	MVector temp(rhs);
		for (int i=0; i<temp.size(); i++) 
		{	if(temp.size() != lhs.size())
			{
				std::cout<<"ERROR: Vectors being summed must be of the same size" <<std::endl;
				abort();
			}
			else temp[i]+=lhs[i];
		}
		return temp;

}

// Operator overload for "vector - vector"
inline MVector operator-(const MVector& lhs, const MVector& rhs)
{
	MVector temp(lhs);
		for (int i=0; i<temp.size(); i++) 
		{	if(temp.size() != rhs.size())
			{
				std::cout<<"ERROR: Vector being subtracted must be of the same size" <<std::endl;
				abort();
			}
			else temp[i]-=rhs[i];
		}
		return temp;

}


//  Operator overload to output a vector to the screen or file
inline std::ostream& operator<<(std::ostream& os, const MVector& v)
{
	int n = v.size();
	os << "(";
	for(int i = 0; i<n; i++)
	{
		if (i != n-1) os << v[i] <<", ";	
		
		else os << v[i];
	}
	os << ")";

	return os;
}