#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <time.h>
#include <sstream> //for osstream
#include <algorithm>

#include "Message.h"
#include "Point.h"

using namespace std;


std::string Point::BackSlash = "/";
std::string Point::DBext = Message::GetDBext();
std::string Point::PointDatabase = Message::GetPointDatabase();
std::string Point::FullResRootName = Message::GetFullResRootName();
std::string Point::CurrentPointDatabase = Message::GetCurrentPointDatabase();
std::string Point::PointFolderRootName = Message::GetPointFolderRootName();
std::string Point::FunResFolderRootName = Message::GetFunResFolderRootName();
std::string Point::FunResRootName = Message::GetFunResRootName();
std::string Point::PointResRootName = Message::GetPointResRootName();


//constructor
Point::Point(int id, double xi, double y){
  m_id = id;
  m_xi = xi;
  m_y = y;
  m_MC.resize(Message::GetNFUN());
  m_MC_to_do.resize(Message::GetNFUN());
  m_NResFiles.resize(Message::GetNFUN());
  m_average.resize(Message::GetNFUN());
  m_stddev.resize(Message::GetNFUN());
  std::stringstream iid;
  iid << m_id;
  m_id_str = iid.str();
  m_myDir = Message::GetResDir() + PointFolderRootName + m_id_str + BackSlash;
}

Point::Point(Point *p){
  m_IdDir = p->GetIdDir();
  m_id = p->GetId();
  m_xi = p->GetXi();
  m_y = p->GetY();
  m_MC.resize(Message::GetNFUN());
  m_MC_to_do.resize(Message::GetNFUN());
  m_NResFiles.resize(Message::GetNFUN());
  m_average.resize(Message::GetNFUN());
  m_stddev.resize(Message::GetNFUN());

  for(int i = 0; i < Message::GetNFUN(); i++)
    {
      m_MC[i]       = p->GetMC(i);
      m_MC_to_do[i] = p->GetMCToDo(i);
      m_NResFiles[i]= p->GetNResFiles(i);
      m_average[i]  = p->GetAverage(i);
      m_stddev[i]   = p->GetStdDev(i);
    }
}

void Point::SetMCToDo(std::vector<int> *desired_MC)
{
  for (int ifun = 0; ifun < Message::GetNFUN(); ifun ++)
    {
      if((*desired_MC)[ifun] - m_MC[ifun] > 0)
	m_MC_to_do[ifun] = (*desired_MC)[ifun] - m_MC[ifun];
      else
	m_MC_to_do[ifun] = 0;
      Message::Debug("SetMCToDo, Point Id %d, ifun = %d, MC done= %d, MC desired = %d, MC to do = %d", m_id, ifun, m_MC[ifun], (*desired_MC)[ifun], m_MC_to_do[ifun]);
    }  
}


void Point::Print()
{
  Message::Info("Printing point informations...");
  Message::Info("id = %d", m_id);
  Message::Info("xi = %g", m_xi);
  Message::Info("y  = %g", m_y);
}


void Point::ReadAllPoints(std::vector<Point*> *PointDone)
{
  //Read all folder and build the right Point
}

void Point::CreatePointsToDo(std::vector<Point*> *PointToDo, std::vector<Point*> *PointDone)
{
  //Use the param file to build point according to :
  // -- The given grid (param)
  // -- The already existing point (PointDone)
  //Destruction of the point must be achieved by the user !
}

void Point::LaunchMC()
{
  int MC_MAX = 0;
  const int NFUN = Message::GetNFUN();
  for (int ifun = 0 ; ifun < NFUN ; ifun ++)
    MC_MAX = std::max(MC_MAX, m_MC_to_do[ifun]);
  Message::Info("[Proc %d] I will do %d MC tests on point with id %d and (xi,y) = (%g, %g)", Message::GetRank(), MC_MAX, m_id, m_xi, m_y);
  //Prepare Res file
  //#pragma omp parallel for private(imc)
  // NFUN vectors of different sizes containing the results...
  std::vector<std::vector<double>* > resultsMC(NFUN);
  for (int ifun = 0; ifun < NFUN; ifun ++)
    {
      resultsMC[ifun] = new std::vector<double>;
      resultsMC[ifun]->reserve(MC_MAX); //Avoiding memory problem
    }
  for (int imc = 0 ; imc < MC_MAX ; imc++)
    {
      std::vector<double> res_int;
      ShortCyclePlus(&res_int);
      for (int ifun = 0; ifun < NFUN; ifun ++)
	resultsMC[ifun]->push_back(res_int[ifun]);
    }
  Message::Debug("MC TERMINATED FOR POINT %d", m_id);
  //Updating files
  WriteOnFile(&resultsMC);
  //Cleaning
  for (int ifun = 0; ifun < NFUN; ifun ++)
      delete resultsMC[ifun];
  Message::Info("[Proc %d] Finnished %d MC tests on point %g %g", Message::GetRank(), MC_MAX, m_xi, m_y);
}

void Point::ShortCyclePlus(std::vector<double> *integrals)
{
  int NFUN = Message::GetNFUN();
  integrals->resize(NFUN);
  for(int i = 0; i < NFUN; i++)
    (*integrals)[i] = 0.0;
  //run a trajectory starting from (xi,y) \in DELTAPLUS of the solution on the boundary.
  int cpt = 0;
  double t = 0.;
  double T = Message::GetT();
  double dt = Message::GetDt();
  double sdt = Message::GetSdt();
  double ro = Message::GetRo();
  double roc = Message::GetRoc();
  double c0 = Message::GetC();
  double k = Message::GetK();
  double Y = Message::GetY();
  double alpha = Message::GetAlpha();
  double lambda = Message::GetLambda();
  int it_max = T/dt;
  double xi = m_xi, y = m_y;
  for(int it = 0 ; it < it_max ; it++)
    {
      t += dt;
      if(t > T){ t = T;}
      //noise
      double g1 = gauss();
      double g2 = gauss();
      double xi_aux = xi, y_aux = y;
      xi += -alpha*xi_aux*dt + sdt*g1;
      y  += -(alpha*xi_aux + c0*y_aux + k*Y)*dt + sdt*(ro*g1 + roc*g2);
      //Compute integrals (Hard coded !)
      for(int i=0; i < NFUN ; i++)
	(*integrals)[i] += dt*exp(-lambda*t)*gplus(xi, y, i);
      if (y<0)
	{
	  for (int i = 0; i < NFUN ; i++)
	    (*integrals)[i] += exp(-lambda*t)*f(xi, i);
	  break;
	}
    }
  return;
}

double Point::gplus(double xi, double y, int ifun)
{
  if( ifun == 0) return f(xi, ifun)*f(y, ifun);
  else return -1.;
}

double Point::f(double xi, int ifun){ 
  if(ifun == 0) return exp(-xi*xi);
  else return -1.;
}


//Write on file and update data about the point...
void Point::WriteOnFile(std::vector<std::vector<double>*> *results)
{
  int NFUN = Message::GetNFUN();
  //Index of the functions that have new results to be stored !
  std::vector<int> FunWithNewResults;
  FunWithNewResults.reserve(NFUN);
  for (int ifun = 0; ifun < NFUN; ifun ++)
    {
      if(m_MC_to_do[ifun] > 0)
	FunWithNewResults.push_back(ifun);
    }
  int NFUNWithNewRes = FunWithNewResults.size();
  //Build a new res file for this Point and every functions
  std::vector<std::ofstream *> fRes(NFUNWithNewRes);
  for (int ifunAux = 0; ifunAux < NFUNWithNewRes; ifunAux++)
    {
      int ifunId = FunWithNewResults[ifunAux];
      Message::Debug("ifunId = %d", ifunId);
      std::ostringstream osifun, osiNbFiles;
      osifun << ifunId;
      osiNbFiles << m_NResFiles[ifunId];
      std::string resFileName = m_myDir + FunResFolderRootName + osifun.str() + BackSlash + PointResRootName + osiNbFiles.str()  + DBext;
      Message::Debug("Writing on %s",resFileName.c_str());
      fRes[ifunAux] = new ofstream(resFileName.c_str(), std::ios_base::out); 
      if(!fRes[ifunAux]->is_open()) Message::Warning("Problem opening file \"%s\"", resFileName.c_str());
    }

  //THIS IS JUST TO DELETE SOME RESULTS TO GET THE EXACT NUMBER OF MC... YES, IT'S WEIRD !!
  std::vector<int> HowManyRes(NFUNWithNewRes);
  for (int ifunAux = 0; ifunAux < NFUNWithNewRes; ifunAux++)
    {
      int ifunId = FunWithNewResults[ifunAux];
      int sizeResults = (*results)[ifunId]->size();
      HowManyRes[ifunAux] = std::max(sizeResults, m_MC_to_do[ifunId]);
      //The first number if the number of results stored in the file
      *(fRes[ifunAux]) << HowManyRes[ifunAux] << "\n";	  
    }

  //Write on files
  for (int ifunAux = 0; ifunAux < NFUNWithNewRes; ifunAux++)
    {
      int ifunId = FunWithNewResults[ifunAux];
      for (int ires =0; ires < HowManyRes[ifunAux]; ires++)
	{
	  *(fRes[ifunAux]) << (*(*results)[ifunId])[ires] << "\n";	  
	}
      fRes[ifunAux]->close();
    }

  //Update data about the points
  for (int ifunAux = 0; ifunAux < NFUNWithNewRes; ifunAux ++)
    {
      int ifunId = FunWithNewResults[ifunAux];
      m_NResFiles[ifunId] ++; 
      m_MC[ifunId] += m_MC_to_do[ifunId];
      m_MC_to_do[ifunId] = 0;
    }

  //Cleaning
  for (int ifunAux = 0; ifunAux < NFUNWithNewRes; ifunAux ++)
    delete fRes[ifunAux];

}
