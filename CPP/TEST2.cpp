#include <vector>
#include <iostream>
#include <math.h>


std::vector<double> tempLeft;
std::vector<double> tempRight;
std::vector<double> tempUp;
std::vector<double> tempDown;
std::vector<double> templength;
std::vector<double> residual;

double norm(std::vector<double>& u, std::vector<double>& v)
{
	double temp = 0;
	for(int i = 0; i< u.size(); i++) temp += (u[i] - v[i])*(u[i] - v[i]);
	double norm = sqrt(temp);
	return norm;
}

class Node
{
public:
	bool fixed;
	double x, y;
	int eqnindx, eqnindy;
	std::vector<int> neighbours;
	std::vector<int> linkNeighbours;
};

class Link
{
public:
	bool vertical;
	int aindex, bindex;
	double length;
};

class Mesh
{
public:
	std::vector<Node> nodes;
	std::vector<Link> links;

	void Create(int nx, int ny)
	{
		int eqnno = 0;

		nodes = std::vector<Node>(nx*ny);

		for (int x=0; x<nx; x++) // assign node coordinates and indices
		{
			for (int y=0; y<ny; y++)
			{
				int thisNode = x+y*nx;
				nodes[thisNode].x = x;
				nodes[thisNode].y = y;
				if (y>0 && y<(ny-1)) // free nodes
				{
					nodes[thisNode].fixed = false;
					nodes[thisNode].eqnindx = eqnno;
					eqnno++;
					nodes[thisNode].eqnindy = eqnno;
					eqnno++;
				}
				else // fixed nodes
				{
					nodes[thisNode].fixed = true;
					nodes[thisNode].eqnindx = -1;
					nodes[thisNode].eqnindy = -1;
				}
			}
		}

		for (int i=0; i<nx*ny; i++) // assign neighbours
		{

			if (nodes[i].x == 0) // left side free nodes
			{
				nodes[i].neighbours.resize(3);
				nodes[i].neighbours[0] = i + 1;
				nodes[i].neighbours[1] = i + nx;
				nodes[i].neighbours[2] = i - nx;
			} 
			else if (nodes[i].x == nx-1)  // right side free ndoes
			{	
				nodes[i].neighbours.resize(3);
				nodes[i].neighbours[0] = i - 1;
				nodes[i].neighbours[1] = i + nx;
				nodes[i].neighbours[2] = i - nx;
			}
			else  // "middle" free nodes with 4 neighbours
			{	
				nodes[i].neighbours.resize(4);
				nodes[i].neighbours[0] = i - 1;
				nodes[i].neighbours[1] = i + 1;
				nodes[i].neighbours[2] = i + nx;
				nodes[i].neighbours[3] = i - nx;
			}
		}

		links = std::vector<Link>(ny*(nx-1) + nx*(ny-1));
		int linkIndex = 0;
		int nLinks = nx*(ny-1) + ny*(nx-1);

		for (int n=0; n<nLinks; n++)
		{
			if (n<nx*(nx-1)) links[n].vertical = false;
			else links[n].vertical = true;
		}

		for (int x=0; x<(nx-1); x++) // Horizontal links
		{
			for (int y=0; y<ny; y++)
			{
				links[linkIndex].aindex = x+y*nx;
				links[linkIndex].bindex = x+1+y*nx;
				links[linkIndex].length = 1;
				linkIndex++;
			}
		}
		for (int x=0; x<nx; x++)
		{
			for (int y=0; y<(ny-1); y++) // Vertical links
			{
				links[linkIndex].aindex = x+y*nx;
				links[linkIndex].bindex = x+(y+1)*nx;
				links[linkIndex].length = 1;
				linkIndex++;
			}
		}


		for (int thisNode=0; thisNode<nx*ny; thisNode++)
		{	
			int linkNeighboursIndex = 0;

			if (nodes[thisNode].fixed == false)
			{
				if (nodes[thisNode].x > 0 && nodes[thisNode].x < nx-1)
				{
					nodes[thisNode].linkNeighbours.resize(4);
				}
				else
				{
					nodes[thisNode].linkNeighbours.resize(3);
				}

				for (int linkIndex=0; linkIndex<nLinks; linkIndex++)
				{
					if (links[linkIndex].aindex == thisNode || links[linkIndex].bindex == thisNode)
					{	
						nodes[thisNode].linkNeighbours[linkNeighboursIndex] = linkIndex;
						linkNeighboursIndex++;
					}
				}
			}
		}
	}
};

int main()
{
	int nx = 4;
	int ny = 4;
	int nLinks = nx*(ny-1) + ny*(nx-1);
	int K = 3; // link constant
	Mesh m;
	m.Create(nx,ny);

	double stretch = 2;   // CONTROL PARAMETER


	// STORE TEMPORARY LENGTH OF LINKS HERE
	templength = std::vector<double>(nLinks);

	for (int linkIndex = 0; linkIndex < nLinks; linkIndex++) 
	{
		templength[linkIndex] = m.links[linkIndex].length;
	}

	
	int y = ny-1; // Top row

	for(int x=0; x<nx; x++) // stretching the top end
	{
		int thisNode = x+y*nx;
		m.nodes[thisNode].y += stretch;
	}

	for (int linkIndex = 0; linkIndex < nLinks; linkIndex++) // calculating new lengths
	{
		if (m.links[linkIndex].vertical == true)
		{
			tempUp = std::vector<double>(2);
			tempDown = std::vector<double>(2);
			tempUp[0] = m.nodes[m.links[linkIndex].bindex].x;
			tempUp[1] = m.nodes[m.links[linkIndex].bindex].y;
			tempDown[0] = m.nodes[m.links[linkIndex].aindex].x;
			tempDown[1] = m.nodes[m.links[linkIndex].aindex].y;
			m.links[linkIndex].length = norm(tempDown, tempUp);
		}
		else
		{
			tempLeft = std::vector<double>(2);
			tempRight = std::vector<double>(2);
			tempLeft[0] = m.nodes[m.links[linkIndex].aindex].x;
			tempLeft[1] = m.nodes[m.links[linkIndex].aindex].y;
			tempRight[0] = m.nodes[m.links[linkIndex].bindex].x;
			tempRight[1] = m.nodes[m.links[linkIndex].bindex].y;
			m.links[linkIndex].length = norm(tempRight, tempLeft);
		}
	}


	// STORE ALL RESIDUALS HERE...... HOW???


	// residual = std::vector<double>(2*(nx*ny - 2*(nx-1)));

	// for ( int eqnno = 0; eqnno<residual.size(); eqnno += 2)
	// {
	// 	for (int thisNode=0; thisNode<nx*ny; thisNode++)
	// 	{
	// 		for (int neighbour = 0; neighbour < m.nodes[thisNode].linkNeighbours.size(); neighbour++)
	// 		{	
	// 			int thisLink = m.nodes[thisNode].linkNeighbours[neighbour];
	// 			Node thisLinkA = m.nodes[m.links[thisLink].aindex];
	// 			Node thisLinkB = m.nodes[m.links[thisLink].bindex];
	// 			if (m.links[thisLink].aindex == thisNode)
	// 			{
	// 				residual[eqnno] += K*(m.links[thisLink].length - templength[thisLink])*((thisLinkB.x-m.nodes[thisNode].x)/(m.links[thisLink].length));
	// 				residual[eqnno+1] += K*(m.links[thisLink].length - templength[thisLink])*((thisLinkB.y-m.nodes[thisNode].y)/(m.links[thisLink].length));
	// 			}
	//	 			else if (m.links[thisLink].bindex == thisNode)
	//	 			{
	//	 				residual[eqnno] += K*(m.links[thisLink].length - templength[thisLink])*((thisLinkA.x-m.nodes[thisNode].x)/(m.links[thisLink].length));
	//	 				residual[eqnno+1] += K*(m.links[thisLink].length - templength[thisLink])*((thisLinkA.y-m.nodes[thisNode].y)/(m.links[thisLink].length));
	//	 			}
	//	 		}
	//	 	}
	//	 }
	//
	//	 std::cout<< residual[20]<<'\n';



}
