#include <iostream>
using namespace std;
#include <math.h>
#include <vector>
#include "myFunctions.h"


/// Forward declaration of Node class as Node pointers are 
/// required in Spring class
class Node; 

//===================================================================
//===================================================================
//===================================================================
/// A spring is connected to two nodes
/// TODO add spriing stiffness
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


 /// access function to the stiffness of the spring
 double& stiffness()
 {
  return K;
 }

 double& strain();


 /// Set the length of the spring
void set_length();

/// Read the length of the spring
double get_length()
{
  return length;
}

/// is the spring vertical?
bool vertical;



private:
 /// Pointers to the nodes at the start and end of the spring. 
 vector<Node*> Node_pt;
 /// Length of the spring
 double length;
 /// stiffness of the spring
 double K;
 /// how much the spring has been stretched
 double s;
};

//===================================================================
//===================================================================
//===================================================================
/// Node has a spatial position (unknown and therefore associated
/// equation number for each of its coordinates) and is connected
/// to a certain number of springs.
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

double& Spring::strain()
{
  return s;
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

//==========================================================
/// Create the mesh using the Spring and Node classes
//==========================================================

class Mesh
{
 public:
  /// constructor (nX contains number of nodes in each direction)
  Mesh(const vector<unsigned>& nX)
    : nodes(num_nodes(nX), nX.size()), springs(num_springs(nX))
    { 
      assign_coordinates(nodes, nX); 
      assign_neighbours(springs, nX);
      dim = nX.size();
      nx = nX;
    }

	/// Access function to the ith node of the mesh.
	Node& node(const unsigned& i)
	{
		return nodes[i];
	}

  /// Access function to the ith spring of the mesh.
  Spring& spring(const unsigned& i)
  {
    return springs[i];
  }

  /// function to set and get the number of equations
  unsigned& num_equations();

  /// function to displace the nodes at the top of the mesh
  void stretch(double e);

  /// function to set the spring stiffnesses
  void stiffness(double primary, double secondary);

  /// Function declaration for assigning coordinates to nodes and
  /// equation numbers to coordinates
  void assign_coordinates(std::vector<Node>& nodes, std::vector<unsigned> nX);
  /// Function declaration for assigning end nodes to the springs
  /// and also to tell the nodes which springs are attached to them.
  void assign_neighbours(std::vector<Spring>& springs, std::vector<unsigned> nX);
  /// Function to get the residuals to zero for all the "free" nodes
  /// e: stretch
  void get_residuals(double e);
  /// Function to get the residuals with perturbed coordinates ( for the Jacobian )
  /// e: stretch, eps: perturbation
  void get_residuals_perturbed(double e, double eps);




private:
  /// Declare vectors to hold the nodes and springs.
  vector<Node> nodes;
  vector<Spring> springs;
  unsigned dim;
  vector<unsigned> nx;
  /// Declare vector to store the residuals
  //vector<double> res;
  unsigned num_of_equations;
};

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

// Member functions

//==============================================================
/// Set the pointer to the i-th node and tell the node that it's
/// connected to the current spring!
//==============================================================

unsigned& Mesh::num_equations()
{
  return num_of_equations;
}

void Spring::set_node_pt(const unsigned& i, Node* node_pt)
{
 Node_pt[i]=node_pt;
 Node_pt[i]->add_spring_pt(this);
}

//==============================================================
/// Function to calculate the length of a spring. "diff_norm" defined
/// in "myFunctions.cpp".
//==============================================================
void Spring::set_length()
{
  length = diff_norm(Node_pt[0]->get_vector(), Node_pt[1]->get_vector());
}

//==============================================================
/// Assign coordinates to the nodes. Input a vector of nodes and 
/// the spatial dimensions.
//==============================================================
void Mesh::assign_coordinates(std::vector<Node>& nodes, std::vector<unsigned> nX)
{
  unsigned dim = nX.size();

  unsigned nx = nX[0];
  unsigned ny = nX[1];
  unsigned nz = nX[2];
  // node index
  unsigned thisNode;
  // 2d mesh
  if (dim == 2)
  { 
    int eqnno = 0;
    // loop through y coordinate
    for (unsigned y=0; y<ny; y++)
    {
      // loop through x coordinate
      for (unsigned x=0; x<nx; x++)
      {
        // update node index
        thisNode = x + y*nx;    
        // member access functions to node coordinates  
        nodes[thisNode].x(0) = x;
        nodes[thisNode].x(1) = y;
        // Assign equation numbers to free nodes
        if (y > 0 && y < ny-1)
        {
          nodes[thisNode].eqn_number(0) = eqnno;
          eqnno++;
          nodes[thisNode].eqn_number(1) = eqnno;
          eqnno++;  
        }
        // Pin the fixed nodes
        else
        {
          nodes[thisNode].pin(0);
          nodes[thisNode].pin(1);
        }
      }
    }
    num_equations() = eqnno;
    std::cout << "The number of equations is " << num_equations() << std::endl;
  }

  // 3d mesh
  else if (dim == 3)
  { 
    int eqnno = 0;
    // loop through z coordinate
    for (unsigned z=0; z<nz; z++)
    {
      // loop through y coordinate
      for (unsigned y=0; y<ny; y++)
      {
        // loop through x coordinate
        for (unsigned x=0; x<nx; x++)
        {
          // update node index
          thisNode = x+y*nx+z*(ny*nx);
          // member access functions to node coordinates
          nodes[thisNode].x(0) = x;
          nodes[thisNode].x(1) = y;
          nodes[thisNode].x(2) = z;
          /// assign equation numbers to the free nodes
          if (z > 0 && z < nz-1)
          {
            nodes[thisNode].eqn_number(0) = eqnno;
            eqnno++;
            nodes[thisNode].eqn_number(1) = eqnno;
            eqnno++;
            nodes[thisNode].eqn_number(2) = eqnno;
            eqnno++;
          }
          /// pin the fixed nodes
          else
          {
            nodes[thisNode].pin(0);
            nodes[thisNode].pin(1);
            nodes[thisNode].pin(2);
          }
        }
      }
    }
    num_equations() = eqnno;
  }
}

  //==============================================================
  /// Tell each spring what two nodes it is connected to and also
  /// tell each node what springs are attached to it.
  //==============================================================
  void Mesh::assign_neighbours(std::vector<Spring>& springs, std::vector<unsigned> nX)
  {
    unsigned dim = nX.size();
    unsigned nx = nX[0];
    unsigned ny = nX[1];
    unsigned nz = nX[2];
    /// spring index
    unsigned thisSpring;
    /// node index
    unsigned thisNode;
    if (dim == 2)
    {
        thisSpring = 0;
        /// springs in the x-direction
        for (unsigned y=0; y<ny; y++)
        {
          for (unsigned x=0; x<nx-1; x++)
          {
            thisNode = x + y*nx;
            springs[thisSpring].set_node_pt(0, &nodes[thisNode]);
            springs[thisSpring].set_node_pt(1, &nodes[thisNode+1]);
            springs[thisSpring].set_length();
            springs[thisSpring].vertical = false;
            thisSpring++;
          }
        }
        /// springs in the y-direction
        for (unsigned y=0; y<ny-1; y++)
        {
          for (unsigned x=0; x<nx; x++)
          {
            thisNode = x + y*nx;
            springs[thisSpring].set_node_pt(0, &nodes[thisNode]);
            springs[thisSpring].set_node_pt(1, &nodes[thisNode + nx]);
            springs[thisSpring].set_length();
            springs[thisSpring].vertical = true;
            thisSpring++;
          }
        }
    }
    else if (dim == 3)
    {
      thisSpring = 0;
      /// springs in the x-direction
      for (unsigned z=0; z<nz; z++)
      {
        for(unsigned y=0; y<ny; y++)
        {
          for (unsigned x=0; x<nx-1; x++)
          {
            thisNode = x+y*nx+z*(ny*nx);
            springs[thisSpring].set_node_pt(0, &nodes[thisNode]);
            springs[thisSpring].set_node_pt(1, &nodes[thisNode+1]);
            springs[thisSpring].set_length();
            springs[thisSpring].vertical = false;
            thisSpring++;

          }
        }
      }
      /// springs in the y-direction
      for (unsigned z=0; z<nz; z++)
      {
        for (unsigned y=0; y<ny-1; y++)
        {
          for (unsigned x=0; x<nx; x++)
          {
            thisNode = x+y*nx+z*(ny*nx);
            springs[thisSpring].set_node_pt(0, &nodes[thisNode]);
            springs[thisSpring].set_node_pt(1, &nodes[thisNode + nx]);
            springs[thisSpring].set_length();
            springs[thisSpring].vertical = false;
            thisSpring++;
          }
        }
      }
      /// springs in the z-direction
      for (unsigned z=0; z<nz-1; z++)
      {
        for (unsigned y=0; y<ny; y++)
        {
          for (unsigned x=0; x<nx; x++)
          {
            thisNode = x+y*nx+z*(ny*nx);
            springs[thisSpring].set_node_pt(0, &nodes[thisNode]);
            springs[thisSpring].set_node_pt(1, &nodes[thisNode + (nx*ny)]);
            springs[thisSpring].set_length();
            springs[thisSpring].vertical = true;
            thisSpring++;
          }
        }
      }
    }
  }

  //==============================================
  /// function to stretch the mesh at one end by e
  //==============================================
  void Mesh::stretch(double e)
  {
    if (nx.size() == 2)
    {
      /// only stretch the top nodes
      unsigned y = nx[1]-1;
      /// loop over top row
      for (unsigned x=0; x<nx[0]; x++)
      {
        unsigned thisNode = x+y*nx[0];
        nodes[thisNode].x(1) += e;  
      }

    }
    else if (nx.size() == 3)
    {
      unsigned z = nx[2] - 1;

      for (unsigned y=0; y<nx[1]; y++)
      {
        for (unsigned x=0; x<nx[0]; x++)
        {
          unsigned thisNode = x+y*nx[0]+z*(nx[1]*nx[0]);
          nodes[thisNode].x(2) += e;
        }
      }
    }
  }

  //==============================================================
  /// Function to set the residuals of each node in the mesh.
  //==============================================================
  void Mesh::get_residuals(double e)
  { 
    /// stretch the mesh
    stretch(e);
    /// vector to hold the residuals.
    vector<double> res(num_equations());
    /// vector to hold the perturbed residuals
    vector<double> res_perturbed(num_equations()*dim);
    
    /// loop over springs
    for (int i=0; i<springs.size(); i++)
    {
      /// this spring
      Spring* thisSpring = &springs[i];
      /// store resting length
      double oldLength = thisSpring->get_length();
      /// calculate new length 
      thisSpring->set_length();
      /// store new length
      double newLength = thisSpring->get_length();
      std::cout << "old length is " << oldLength << " and new length is " << newLength << std::endl;

      /// Check if the nodes at the ends of each spring are pinned.
      /// If they're not pinned, add the force from the spring the 
      /// respective residual equation.
      if (thisSpring->node_pt(0)->is_pinned(0) == false)
      {
        res[thisSpring->node_pt(0)->eqn_number(0)] += thisSpring->stiffness()*((newLength - oldLength)*((thisSpring->node_pt(1)->x(0) - thisSpring->node_pt(0)->x(0))))/(newLength);
      }
      if (thisSpring->node_pt(0)->is_pinned(1) == false)
      {
        res[thisSpring->node_pt(0)->eqn_number(1)] += thisSpring->stiffness()*((newLength - oldLength)*((thisSpring->node_pt(1)->x(1) - thisSpring->node_pt(0)->x(1))))/(newLength);
      }
      if (thisSpring->node_pt(1)->is_pinned(0) == false)
      {
        res[thisSpring->node_pt(1)->eqn_number(0)] -= thisSpring->stiffness()*((newLength - oldLength)*((thisSpring->node_pt(1)->x(0) - thisSpring->node_pt(0)->x(0))))/(newLength);
      }
      if (thisSpring->node_pt(1)->is_pinned(1) == false)
      {
         res[thisSpring->node_pt(1)->eqn_number(1)] -= thisSpring->stiffness()*((newLength - oldLength)*((thisSpring->node_pt(1)->x(1) - thisSpring->node_pt(0)->x(1))))/(newLength);              
      }
    }
    /// Test output
    for (int i=0; i<res.size(); i++)
    {
      std::cout<< "the residual is "<<res[i]<<std::endl;
    }
  }


    /// function to set the spring stiffnesses
  void Mesh::stiffness(double primary, double secondary)
  {
    for (int i=0; i<springs.size(); i++)
    {    
        /// check if spring is vertical or horizontal
        if (springs[i].vertical == true)
        {
          springs[i].stiffness() = primary;
        }
        else
        {
          springs[i].stiffness() = secondary;
        }
    }
  }


///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

int main()
{
  /// spatial dimensions
  unsigned nx = 3;
  unsigned ny = 3;
  /// store into vector
  vector<unsigned> nX(2);
  nX[0] = nx;
  nX[1] = ny;
  /// create a mesh
  Mesh m(nX);
  /// set the spring stiffnesses (primary, secondary)
  m.stiffness(3, 2);
  m.get_residuals(5.0);
};
