#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <math.h>
#include <vector>
#include <algorithm>

#include "Message.h"
#include "Mesh.h"

Mesh::Mesh(std::vector<double> X,  std::vector<double> Y)
{
  int npx = X.size();
  int npy = Y.size();
  if(npx != npy)
    {
      Message::Warning("X and Y not of the same size, Mesh impossible, abording...");
      Message::Finalize(EXIT_FAILURE);
    }
  m_X.resize(npx);
  m_Y.resize(npx);
  m_Z.resize(npx);
  for(int i = 0; i< npx; i++)
    {
      m_X[i] = X[i];
      m_Y[i] = Y[i];
      m_Z[i] = 0;
    }
  //juste in case of ...
  m_connectivity.reserve(npx*4);
  /*  m_res.resize(Message::GetNFUN());
  for (int i =0; i<Message::GetNFUN();i++)
  m_res[i].resize(npx);*/
}

void Mesh::SetRes(std::vector<std::vector<double> > *res)
{
  /*  int npx = m_X.size();
  for (int ifun =0; i<Message::GetNFUN();i++)
      for(int i = 0; i< npx; i++)
      m_res[ifun][i] = res[ifun][i];*/
  m_res = res;
}


void Mesh::Update()
{
  ReArrange();

  m_connectivity.resize(2);
  m_connectivity[0].resize(3);
  m_connectivity[0][0] = 0;
  m_connectivity[0][1] = 1;
  m_connectivity[0][2] = 2;
  m_connectivity[1].resize(3);
  m_connectivity[1][0] = 0;
  m_connectivity[1][1] = 2;
  m_connectivity[1][2] = 3;
  int Npoints = m_X.size();
  for(int i = 4; i < Npoints; i++)
    {
      std::vector<int> res;
      IsOnAnEdge(i, &res);
      int sizeres = res.size();
      if(sizeres>0)
	{
	  for (int j = 0 ; j < sizeres; j+=2)
	    {
	      SplitTriangleFromEdge(res[j], res[j+1], i);	    }
	  continue;
	}
      int iTri = IsInATriangle(i);
      if(iTri > -1)
	{
	  SplitTriangle(iTri, i);
	  continue;
	}
    }
}


void Mesh::ReArrange()
{
  std::vector<int> ind_corner, ind_interior, ind_border;
  SeparatePoints(&ind_corner, &ind_interior, &ind_border);
  int ncorner = ind_corner.size(), ninterior = ind_interior.size(),nborder = ind_border.size();
  //Auxiliary table
  std::vector<double> XX, YY;
  int npoints = m_X.size();
  XX.resize(npoints);
  YY.resize(npoints);
  for(int i =0; i< npoints; i++)
    {
      XX[i] = m_X[i];
      YY[i] = m_Y[i];
      m_X[i] = 0;
      m_Y[i] = 0;
    }
  int cpt = 0;
  if(ncorner !=4)
    {
    Message::Warning("Huge problem, no 4 corners found !! Abording Mesh computation...");
    return;
    }
  for (int i =0; i < ncorner; i++)
    {
      m_X[cpt] = XX[ind_corner[i]];
      m_Y[cpt] = YY[ind_corner[i]];
      cpt++;
    }
  for (int i =0; i < ninterior; i++)
    {
      m_X[cpt] = XX[ind_interior[i]];
      m_Y[cpt] = YY[ind_interior[i]];
      cpt++;
    }
  for (int i =0; i < nborder; i++)
    {
      m_X[cpt] = XX[ind_border[i]];
      m_Y[cpt] = YY[ind_border[i]];
      cpt++;
    }
}

void Mesh::SeparatePoints(std::vector<int> *ind_corner, std::vector<int> *ind_interior, std::vector<int> *ind_border)
{
  int np = m_X.size();
  ind_corner->resize(4);
  ind_interior->reserve(np);
  ind_border->reserve(np);
  double xmin = *std::min_element(m_X.begin(), m_X.end());
  double xmax = *std::max_element(m_X.begin(), m_X.end());
  double ymin = *std::min_element(m_Y.begin(), m_Y.end());
  double ymax = *std::max_element(m_Y.begin(), m_Y.end());
  std::vector<int> ind_xmin, ind_xmax, ind_ymin, ind_ymax;
  ind_xmin.reserve(np);
  ind_xmax.reserve(np);
  ind_ymin.reserve(np);
  ind_ymax.reserve(np);
  for (int i = 0 ; i<np; i++)
    {
      if(m_X[i] == xmin)
	ind_xmin.push_back(i);
      if(m_X[i] == xmax)
	ind_xmax.push_back(i);
      if(m_Y[i] == xmin)
	ind_ymin.push_back(i);
      if(m_Y[i] == ymax)
	ind_ymax.push_back(i);
    }
  for (int i = 0; i < np; i++)
    {
      if(m_X[i] == xmin && m_Y[i] == ymin)
	(*ind_corner)[0] = i;
      else if(m_X[i] == xmax && m_Y[i] == ymin)
	(*ind_corner)[1] = i;
      else if(m_X[i] == xmax && m_Y[i] == ymax)
	(*ind_corner)[2] = i;
      else if(m_X[i] == xmin && m_Y[i] == ymax)
	(*ind_corner)[3] = i;
      else if(m_X[i] == xmin || m_X[i] == xmax || m_Y[i] == ymin || m_Y[i] == ymax)
	ind_border->push_back(i);
      else
	ind_interior->push_back(i);
    }
}


int Mesh::IsInATriangle(int ind)
{
  //Check if a point is in a triangle or on the edge of.
  double x = m_X[ind];
  double y = m_Y[ind];
  int nTri = m_connectivity.size();
  for (int iT = 0; iT < nTri; iT ++)
    {
      double x0,x1,x2,y0,y1,y2,y3;
      x0 = m_X[m_connectivity[iT][0]]; y0 = m_Y[m_connectivity[iT][0]];
      x1 = m_X[m_connectivity[iT][1]]; y1 = m_Y[m_connectivity[iT][1]];
      x2 = m_X[m_connectivity[iT][2]]; y2 = m_Y[m_connectivity[iT][2]];
      double Area = abs(0.5*(-y1*x2 + y0*(-x1 + x2) + x0*(y1 - y2) + x1*y2));
      double s = 1/(2*Area)*(y0*x2 - x0*y2 + (y2 - y0)*x + (x0 - x2)*y);
      double t = 1/(2*Area)*(x0*y1 - y0*x1 + (y0 - y1)*x + (x1 - x0)*y);
      if(s > 0 && t > 0 && (1-s-t)>0)
	return iT;
    }
  return -1;
}

void Mesh::IsOnAnEdge(int ind, std::vector<int > *res)
{
  //Check if a point is in a triangle or on the edge of one or two triangles
  double x = m_X[ind];
  double y = m_Y[ind];
  int npoints = m_X.size();
  int nTri = m_connectivity.size();
  res->reserve(4);
  for(int iT = 0; iT < nTri; iT++)
    {
      double x0,x1,x2,y0,y1,y2,y3;
      x0 = m_X[m_connectivity[iT][0]]; y0 = m_Y[m_connectivity[iT][0]];
      x1 = m_X[m_connectivity[iT][1]]; y1 = m_Y[m_connectivity[iT][1]];
      x2 = m_X[m_connectivity[iT][2]]; y2 = m_Y[m_connectivity[iT][2]];
      if(IsAligned(x,y,x0,y0,x1,y1))
	{
	  res->push_back(iT);
	  res->push_back(2);
	  continue;
	}
      else if(IsAligned(x,y,x1,y1,x2,y2))
	{
	  res->push_back(iT);
	  res->push_back(0);
	  continue;
	}
      else if(IsAligned(x,y,x0,y0,x2,y2))
	{
	  res->push_back(iT);
	  res->push_back(1);
	  continue;
	}
    }
}

bool Mesh::IsAligned(double x0, double y0, double x1, double y1, double x2, double y2)
{
  double Area = abs(0.5*(-y1*x2 + y0*(-x1 + x2) + x0*(y1 - y2) + x1*y2));
  double scalar_product = (x1-x0)*(x2-x0) + (y1-y0)*(y2-y0);
  return (Area==0 && scalar_product <=0);
}

bool Mesh::IsAVertex(int ind)
{
  int nT = m_connectivity.size();
  for (int i = 0; i < nT; i++)
    {
      if(m_connectivity[i][0]==ind ||m_connectivity[i][1]==ind ||m_connectivity[i][2]==ind)
	return true;
    }
  return false;
}

void Mesh::SplitTriangle(int iTriangle, int iPoint)
{
  int a,b,c;
  a = m_connectivity[iTriangle][0];
  b = m_connectivity[iTriangle][1];
  c = m_connectivity[iTriangle][2];
  //Modify current triangle
  m_connectivity[iTriangle][0] = a;
  m_connectivity[iTriangle][1] = b;
  m_connectivity[iTriangle][2] = iPoint;
  //Add Two new triangles
  int sizeConn = m_connectivity.size();
  m_connectivity.resize(sizeConn+2);
  m_connectivity[sizeConn].resize(3);
  m_connectivity[sizeConn+1].resize(3);
  m_connectivity[sizeConn][0] = b;
  m_connectivity[sizeConn][1] = c;
  m_connectivity[sizeConn][2] = iPoint;
  m_connectivity[sizeConn+1][0] = c;
  m_connectivity[sizeConn+1][1] = a;
  m_connectivity[sizeConn+1][2] = iPoint;
}


void Mesh::SplitTriangleFromEdge(int iTriangle, int Vertex, int iPoint)
{
  int connSize = m_connectivity.size();
  m_connectivity.resize(connSize + 1);
  int a,b,c;
  a = m_connectivity[iTriangle][Vertex];
  b = m_connectivity[iTriangle][(Vertex+1)%3];
  c = m_connectivity[iTriangle][(Vertex+2)%3];
  m_connectivity[iTriangle][0] = a;
  m_connectivity[iTriangle][1] = iPoint;
  m_connectivity[iTriangle][2] = b;
  m_connectivity[connSize].resize(3);
  m_connectivity[connSize][0] = iPoint;
  m_connectivity[connSize][1] = a;
  m_connectivity[connSize][2] = c;
}



void Mesh::PrintMesh(std::string Filename)
{
  int Npoints = m_X.size();
  std::ofstream file((Filename +".msh").c_str());
  file << "$MeshFormat"<< std::endl;
  file << "2.2 0 8"<< std::endl;
  file << "$EndMeshFormat"<< std::endl;
  file << "$Nodes"<< std::endl;
  file << Npoints<< std::endl;
  for (int i = 0; i < Npoints; i++)
    file << i+1 << " " <<m_X[i] << " " << m_Y[i]<< " 0" << std::endl;
  file << "$EndNodes"<< std::endl;
  file << "$Elements"<< std::endl;
  file << m_connectivity.size()<< std::endl;
  for(int i = 0; i < m_connectivity.size(); i++)
    file << i+1 << " 2 2 1 1 " << m_connectivity[i][0]+1 << " " << m_connectivity[i][1]+1<< " " << m_connectivity[i][2]+1<< std::endl;
  file << "$EndElements"<< std::endl;
  file.close();
}


void Mesh::PrintRes(std::string RootFileName)
{
  int iTri = m_connectivity.size();
  for (int ifun = 0; ifun < Message::GetNFUN(); ifun++)
    {
      std::stringstream iifun;
      iifun << ifun;
      std::string FileName = RootFileName + iifun.str() + ".pos";
      std::ofstream posFile(FileName.c_str(), std::ios_base::out);
      posFile << "View \"res\" {"<< "\n";

      for (int iT = 0; iT < iTri; iT++)
	{
	  int a = m_connectivity[iT][0];
	  int b = m_connectivity[iT][1];
	  int c = m_connectivity[iT][2];
	  posFile << "ST(" << m_X[a] <<","<< m_Y[a] << "," << m_Z[a];
	  posFile << "," << m_X[b] <<","<< m_Y[b] << "," << m_Z[b];
	  posFile << "," << m_X[c] <<","<< m_Y[c] << "," << m_Z[c];
	  posFile << "){";
	  posFile << (*m_res)[ifun][a]  << ","<<(*m_res)[ifun][b]<< ","<<(*m_res)[ifun][c] << ",0,0,0};\n";
	}
      posFile << "};\n";
      posFile.close();
    }
}
