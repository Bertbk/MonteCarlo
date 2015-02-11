#ifndef _Mesh_H_
#define _Mesh_H_

#include <vector>
#include <string>


/*
This class builds, from two vectors of X and Y coordinates, a triangular mesh joining these point.
The mesh is clearly not of good quality and 
- The points must be in a rectangular domain with four corners (if not, it's not working)
- To use this class, instanciate a Mesh, then run 
-- Update() to build the mesh
-- SetRes(...) to specify the results
-- PrintRes(Filename) to create a .pos file (GetDP format) setting the result on the points.


*/
class Mesh{
 private:
  std::vector<double> m_X;
  std::vector<double> m_Y;
  std::vector<double> m_Z;
  std::vector<std::vector <int> > m_connectivity;
  std::vector<std::vector<double> > *m_res; // Results on every point (for each functions)

  //Separate corners, border and interior points.
  void SeparatePoints(std::vector<int> *ind_corner, std::vector<int> *ind_interior, std::vector<int> *ind_border);
  //Check if the point of index ind is inside a triangle...
  int IsInATriangle(int ind);
  //or on an (or more) Edge(s) (res[cpt] = index of triangle, res[cpt+1] = which vertex is not on the edge) 
  void IsOnAnEdge(int ind, std::vector<int> *res);
  //Is the point (x,y) on the segment (x0,y0) - (x1,y1) ? 
  bool IsAligned(double x, double y, double x0, double y0, double x1, double y1);
  //Split a triangle in three, where iPoint is a point inside the triangle 
  void SplitTriangle(int iTriangle, int iPoint);
  //Split a triangle in two, where iPoint is on an edge, and Vertex is the only point of the triangle to not be on the edge of iPoint
  //(the two news triangle will have iPoint and Vertex as vertexes (note Vertex =0,1 or 2 (local numbering)))
  void SplitTriangleFromEdge(int iTriangles, int Vertex, int iPoint);

 public:
  Mesh(std::vector<double> X,  std::vector<double> Y);
  void SetRes(std::vector<std::vector<double> > *res);
  //Update compute the 2D triangulation that links all the points, without building new ones.
  // The mesh is probably of bad quality but we do not care: it's only for post processing with GMSH
  void Update();
  //PrintGMSH print the mesh on a file FileName with GMSH format (.msh)
  //Usefull for checking, for example
  void PrintMesh(std::string FileName);
  //Print GetDP .pos format m_res vector on the points.
  void PrintRes(std::string RootFileName);
};


#endif
