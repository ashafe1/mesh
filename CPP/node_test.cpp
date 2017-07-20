using namespace std;
#include <math.h>
#include <iostream>
#include <vector>
#include "myFunctions.h"

class Node2d;
class Node3d;


class Spring2d
{

public:
 /// Constructor
 Spring2d()
  {
   // Each spring is only ever connected to two nodes
   // Also initialise pointers to null.
   Node_pt.resize(2,0);
  }

 /// Read-only function to ith node (i=0 first node and i=1 second node)
 Node2d* node_pt(const unsigned& i) const
  {
   return Node_pt[i];
  }


 /// Set the pointer the i-th node and tell the node that it's
 /// connected to the current spring!
 void set_node_pt(const unsigned& i, Node2d* node_pt);

 /// Set the length of the spring
void set_length();

/// Read the length of the spring
double get_length()
{
  return length;
}


private:
 /// Pointers to the nodes at the start and end of the spring. 
 vector<Node2d*> Node_pt;
 /// Length of the spring
 double length;
};



class Spring3d: public Spring2d
{

public:


 /// Constructor
 Spring3d()
  {
   // Each spring is only ever connected to two nodes
   // Also initialise pointers to null.
   Node_pt.resize(3,0);
  }

 /// Read-only function to ith node (i=0 first node and i=1 second node)
 Node3d* node_pt(const unsigned& i) const
  {
   return Node_pt[i];
  }


 /// Set the pointer the i-th node and tell the node that it's
 /// connected to the current spring!
 void set_node_pt(const unsigned& i, Node3d* node_pt);


private:
 /// Pointers to the nodes at the start and end of the spring. 
 vector<Node3d*> Node_pt;
 /// Length of the spring
 double length;
};

class Node2d
{
public:

 /// (Non-)magic number indicating that the coordinate has not
 /// been classified as pinned or free yet
 static int Not_classified_yet;

 /// (Non-)magic number indicating that the coordinate is pinned
 static int Is_pinned;

 /// Constructor: Pass the spatial dimension
 Node2d()
  {
   // Resize
   X.resize(2,1.0);
   Eqn_number.resize(2,Not_classified_yet);
  }


 /// Function to add a spring to the node
 void add_spring_pt(Spring2d* spring_pt)
  {
   Spring_pt.push_back(spring_pt);
  }

 /// How many springs are attached to this node?
 unsigned nspring()
  {
   return Spring_pt.size();
  }

 /// Access function to the ith spring connected to the node
 Spring2d*& spring_pt(const unsigned& i)
  {
   return Spring_pt[i];
  }
 /// Access function to the position vector of the node
  vector<double>& get_vector()
  {
    return X;
  }

 /// Access function to the coordinates of the node
 double& x(int i)
  {
   return X[i];
  }

 /// Access function to the equation number for each coordinate 
 /// Can be negative if node is pinned in that direction.
 int& eqn_number(const unsigned& i)
  {
   return Eqn_number[i];
  }


 /// Pin the i-th coordinate
 void pin(const unsigned& i)
  {
   Eqn_number[i]=Is_pinned;
  }

 /// Is the i-th coordinate pinned?
 bool is_pinned(const unsigned& i)
  {
   return (Eqn_number[i]==Is_pinned);
  }


private:

 /// Pointers to the springs attatched to the node. 
 vector<Spring2d*> Spring_pt;

 /// Coordinates of the node
 vector<double> X;

 /// Vector containing equation indices for each coordinate direction.
 /// Can be negative if node is pinned in that direction.
 vector<int> Eqn_number;

};

class Node3d: public Node2d
{
public:
	Node3d()
	{
		Node2d();
		X.resize(3,0.0);
	}

 /// Function to add a spring to the node
 void add_spring_pt(Spring3d* spring_pt)
  {
   Spring_pt.push_back(spring_pt);
  }

 /// Access function to the ith spring connected to the node
 Spring3d*& spring_pt(const unsigned& i)
  {
   return Spring_pt[i];
  }
private:

 /// Pointers to the springs attatched to the node. 
 vector<Spring3d*> Spring_pt;

 /// Coordinates of the node
 vector<double> X;

 /// Vector containing equation indices for each coordinate direction.
 /// Can be negative if node is pinned in that direction.
 vector<int> Eqn_number;
};


 /// (Non-)magic number indicating that the coordinate has not
 /// been classified as pinned or free yet
int Node2d::Not_classified_yet=-10;

 /// (Non-)magic number indicating that the coordinate is pinned
int Node2d::Is_pinned=-1;


class Mesh2d
{
 public:
    /// constructor (nX contains number of nodes in each direction)
  Mesh2d(unsigned nx, unsigned ny)
  { 
    /// Function "num_nodes" defined in "myFunctions.cpp" to find the 
    /// total number of nodes.
    unsigned nNodes = nx*ny;
    /// Check the dimension of the mesh and and construct a vector
    /// of the nodes.
    vector<Node2d> nodes(nNodes);
    /// Function "num_springs" defined in "myFunctions.cpp" to find the
    /// total number of springs.
    unsigned nsprings = num_springs(nx, ny);
    /// Vector to hold the springs.
    vector<Spring2d> springs(nsprings);
    /// Function to assign coordinates to all the nodes.
    assign_coordinates(nodes,nx, ny);
    std::cout << nodes[2].x(0) << std::endl;
    std::cout << "dimension of node is "<< nodes[1].get_vector().size()<<std::endl;
  }
	/// Access function to the ith node of the mesh.
	Node2d node(const unsigned& i)
	{
		return nodes[i];
	}

  /// Function declaration for assigning coordinates to nodes
  void assign_coordinates(std::vector<Node2d> nodes, unsigned nx, unsigned ny);

	/// Access function to the ith spring of the mesh.
	Spring2d spring(const unsigned& i)
	{
		return springs[i];
	}

private:
  /// Declare vectors to hold the nodes and springs.
  vector<Node2d> nodes;
  vector<Spring2d> springs;
};

class Mesh3d
{
public:
	Mesh3d(unsigned nx, unsigned ny, unsigned nz)
	{
		/// Function "num_nodes" defined in "myFunctions.cpp" to find the 
	    /// total number of nodes.
	    unsigned nNodes = nx*ny*nz;
	    /// Check the dimension of the mesh and and construct a vector
	    /// of the nodes.
	    vector<Node3d> nodes(nNodes, Node3d());
	    /// Function "num_springs" defined in "myFunctions.cpp" to find the
	    /// total number of springs.
	    unsigned nsprings = num_springs(nx, ny, nz);
	    /// Vector to hold the springs.
	    vector<Spring3d> springs(nsprings);
	    /// Function to assign coordinates to all the nodes.
	    assign_coordinates(nodes,nx, ny, nz);
	    std::cout << nodes[0].x(3) << std::endl;
	    std::cout << "dimension of node is "<< nodes[1].get_vector().size()<<std::endl;
	}
	/// Access function to the ith node of the mesh.
	Node3d node(const unsigned& i)
	{
		return nodes[i];
	}

  /// Function declaration for assigning coordinates to nodes
  void assign_coordinates(std::vector<Node3d> nodes, unsigned nx, unsigned ny, unsigned nz);

	/// Access function to the ith spring of the mesh.
	Spring3d spring(const unsigned& i)
	{
		return springs[i];
	}

private:
  /// Declare vectors to hold the nodes and springs.
  vector<Node3d> nodes;
  vector<Spring3d> springs;
};

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

// Member functions

//==============================================================
/// Set the pointer to the i-th node and tell the node that it's
/// connected to the current spring!
//==============================================================
void Spring2d::set_node_pt(const unsigned& i, Node2d* node_pt)
{
 Node_pt[i]=node_pt;
 Node_pt[i]->add_spring_pt(this);
}

/// Function to calculate the length of a spring. "diff_norm" defined
/// in "myFunctions.cpp".
void Spring2d::set_length()
{
  length = diff_norm(Node_pt[0]->get_vector(), Node_pt[1]->get_vector());
}

/// Assign coordinates to the nodes. Input a vector of nodes and the spatial dimensions
void Mesh2d::assign_coordinates(std::vector<Node2d> nodes, unsigned nx, unsigned ny)
{
    for (unsigned y=0; y<ny; y++)
    {
      for (unsigned x=0; x<nx; x++)
      {
        unsigned thisNode = x + y*nx;      
        nodes[thisNode].x(0) = x;
        nodes[thisNode].x(1) = y;
      }
    }
}
void Mesh3d::assign_coordinates(std::vector<Node3d> nodes, unsigned nx, unsigned ny, unsigned nz)
{
    for (unsigned z=0; z<nz; z++)
    {
      for (unsigned y=0; y<ny; y++)
      {
        for (unsigned x=0; x<nx; x++)
        {
          unsigned thisNode = x+y*nx+z*(ny*nx);
          nodes[thisNode].x(0) = x;
          nodes[thisNode].x(1) = y;
          nodes[thisNode].x(2) = z;
        }
      }
    }
}


///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

int main()
{
	Node2d n;
	Node3d m;
	Mesh3d k(3,3,3);
}