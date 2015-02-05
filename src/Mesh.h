#ifndef _Mesh_H_
#define _Mesh_H_

#include <vector>
#include <string>

class Mesh{
 private:
  std::vector<double> m_X;
  std::vector<double> m_Y;
  std::vector<double> m_Z;
  std::vector<std::vector <int> > m_connectivity;
  std::vector<std::vector<double> > *m_res; // Results on every point (for each functions)

  void ReArrange();
  void SeparatePoints(std::vector<int> *ind_corner, std::vector<int> *ind_interior, std::vector<int> *ind_border);
  int IsInATriangle(int ind);
  void IsOnAnEdge(int ind, std::vector<int> *res);
  bool IsAligned(double x, double y, double x0, double y0, double x1, double y1);
  bool IsAVertex(int ind);
  void SplitTriangle(int iTriangle, int iPoint);
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
  void PrintRes(std::string RootFileName);
};


#endif
