#ifndef _Mesh_H_
#define _Mesh_H_

#include <vector>
#include <string>

class Mesh{
 private:
  std::vector<double> m_X;
  std::vector<double> m_Y;
  std::vector<std::vector <int> > m_connectivity;

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
  void Update();
  void PrintGMSH(std::string FileName);
};


#endif
