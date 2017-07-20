using namespace std;
#include <math.h>
#include <iostream>
#include <vector>
#include "myFunctions.h"


/// Forward declaration of Node class as Node pointers are 
/// required in Spring class
class Node; 
//===================================================================
//===================================================================
//===================================================================
/// A spring is connected to two nodes
// todo add spriing stiffmess
//===================================================================
//===================================================================
//===================================================================
class Spring 
{

public:


 /// Constructor
 Spring()
  {
   // Each spring is only ever connected to two nodes
   // Also initialise pointers to null.
   Node_pt.resize(2,0);
  }

 /// Read-only function to ith node (i=0 first node and i=1 second node)
 Node* node_pt(const unsigned& i) const
  {
   return Node_pt[i];
  }


 /// Set the pointer the i-th node and tell the node that it's
 /// connected to the current spring!
 void set_node_pt(const unsigned& i, Node* node_pt);

 /// Set the length of the spring
void set_length();

/// Read the length of the spring
double get_length()
{
  return length;
}


private:

 /// Pointers to the nodes at the start and end of the spring. 
 /// This needs to be resized to 2 during the mesh creation.
 vector<Node*> Node_pt;
 double length;


};

//===================================================================
//===================================================================
//===================================================================
//===================================================================
//===================================================================
/// Node has a spatial position (unknown and therefore associated
/// equation number for each of its coordinates) and is connected
/// to a certain number of springs.
//===================================================================
//===================================================================
//===================================================================
//===================================================================
//===================================================================
class Node
{ 

public:

 /// (Non-)magic number indicating that the coordinate has not
 /// been classified as pinned or free yet
 static int Not_classified_yet;

 /// (Non-)magic number indicating that the coordinate is pinned
 static int Is_pinned;

 /// Constructor: Pass the spatial dimension
 Node(const unsigned& dim)
  {
   // Resize
   X.resize(dim,0.0);
   Eqn_number.resize(dim,Not_classified_yet);
  }


 /// Function to add a spring to the node
 void add_spring_pt(Spring* spring_pt)
  {
   Spring_pt.push_back(spring_pt);
  }

 /// How many springs are attached to this node?
 unsigned nspring()
  {
   return Spring_pt.size();
  }

 /// Access function to the ith spring connected to the node
 Spring*& spring_pt(const unsigned& i)
  {
   return Spring_pt[i];
  }

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
 vector<Spring*> Spring_pt;

 /// Coordinates of the node
 vector<double> X;

 /// Vector containing equation indices for each coordinate direction.
 /// Can be negative if node is pinned in that direction.
 vector<int> Eqn_number;

};

 /// (Non-)magic number indicating that the coordinate has not
 /// been classified as pinned or free yet
int Node::Not_classified_yet=-10;

 /// (Non-)magic number indicating that the coordinate is pinned
int Node::Is_pinned=-1;

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

//==========================================================
/// Create the mesh using the Spring and Node classes
//==========================================================

class Mesh
{
 public:
//   /// nX is a vector containing the number of nodes in each direction
  Mesh(const vector<unsigned> nX)
  { 
    /// Function "num_nodes" defined in "myFunctions.cpp" to find the 
    /// total number of nodes.
    unsigned nNodes = num_nodes(nX);
    /// Check the dimension of the mesh and enter "dim" into the
    /// constructor for the nodes.
    unsigned dim = nX.size();
    vector<Node> nodes(nNodes,dim);
    /// Function "num_springs" defined in "myFunctions.cpp" to find the
    /// total number of nodes.
    unsigned nsprings = num_springs(nX);
    /// Vector to hold the springs.
    vector<Spring> springs(nsprings);

    assign(nodes,nX);
  }
	/// Access function to the ith node of the mesh.
	Node node(const unsigned& i)
	{
		return nodes[i];
	}

  /// Assign coordinates to the nodes. Input a vector of nodes and the spatial dimensions
 void assign(std::vector<Node> nodes, std::vector<unsigned> nX);

	/// Access function to the ith spring of the mesh.
	Spring spring(const unsigned& i)
	{
		return springs[i];
	}

private:
  /// Declare vectors to hold the nodes and springs.
  vector<Node> nodes;
  vector<Spring> springs;
};

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

// Member functions

//==============================================================
/// Set the pointer the i-th node and tell the node that it's
/// connected to the current spring!
//==============================================================
void Spring::set_node_pt(const unsigned& i, Node* node_pt)
{
 Node_pt[i]=node_pt;
 Node_pt[i]->add_spring_pt(this);
}

void Spring::set_length()
{
  length = diff_norm(Node_pt[0]->get_vector(), Node_pt[1]->get_vector());
}

/// Assign coordinates to the nodes. Input a vector of nodes and the spatial dimensions
void Mesh::assign(std::vector<Node> nodes, std::vector<unsigned> nX)
{
  unsigned dim = nX.size();
  if (dim == 2)
  {
    for (unsigned y=0; y<nX[1]; y++)
    {
      for (unsigned x=0; x<nX[0]; x++)
      {
        unsigned thisNode = x + y*nX[0];      
        nodes[thisNode].x(0) = x;
        nodes[thisNode].x(1) = y;
      }
    }
  }
  else if (dim == 3)
  {
    for (unsigned z=0; z<nX[2]; z++)
    {
      for (unsigned y=0; y<nX[1]; y++)
      {
        for (unsigned x=0; x<nX[0]; x++)
        {
          unsigned thisNode = x+y*nX[0]+z*(nX[1]*nX[0]);
          nodes[thisNode].x(0) = x;
          nodes[thisNode].x(1) = y;
          nodes[thisNode].x(2) = z;
        }
      }
    }
  }
}

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

int main()
{

 
};
