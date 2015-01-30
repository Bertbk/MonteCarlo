#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "Database.h"
#include "Message.h"
#include "Point.h"

std::string Database::DBext = ".db";
std::string Database::PointDatabase = "Points";
std::string Database::FullResRootName = "ResFun";
std::string Database::CurrentPointDatabase = "currentPoint";
std::string Database::FunResRootName = "fun";
std::string Database::PointResRootName = "point_res_";

//Constructor
Database::Database(std::string resdir){
  m_resDir = resdir;
  NpointsDone = 0;
  NpointsToDo = 0;
}

Database::~Database(){
  for (int i = 0; i < NpointsDone; i++)
    delete PointsDone[i];
  for (int i = 0; i < NpointsToDo; i++)
    delete PointsToDo[i];
}


Point* Database::GetPointToDo(int index){
  //Does this work ?
  return PointsToDo[index];
}

void Database::Init()
{
  Message::Info("Init Database...");
  NResByFun.resize(Message::GetNFUN());
  MeanByFun.resize(Message::GetNFUN());
  StdDevByFun.resize(Message::GetNFUN());
  //Read folder then subfolder, then ...
  ParseRootFiles();

}


void Database::ParseRootFiles(){
  Message::Info("ParseRootFiles...");
  //The root directory should have a file points.db and N_FUN files res_funXX.db
  //points.db: contains information about the points (Id, X, Y)
  std::string PointsDbName = Message::GetResDir() + PointDatabase + DBext;
  std::ifstream PointsDb(PointsDbName.c_str(), std::ios_base::in);
  if(!PointsDb.is_open())
    {
      Message::Warning("No database file found... I hope there is no work there because it's gonna be erazed by new results...");
      NpointsDone = 0;
      return;
    }
  else
    {
      Message::Info("Database file found! Let's rock!");
      //      std::string line;
      //      getline(PointsDb, line);
      PointsDb >> NpointsDone;
      MaxId = NpointsDone;
      PointsDone.resize(NpointsDone);
      for (int i = 0; i < NpointsDone; i++)
	{
	  int Id;
	  double xi, y;
	  PointsDb >> Id >> xi >> y;
	  PointsDone[i] = new Point(Id, xi, y);
	}
      PointsDb.close();
    }
  //Now let's attack the other root files, containing the results obtained with differents functions
  for (int ifun = 0; ifun < Message::GetNFUN(); ifun++)
    {
      NResByFun[ifun].resize(MaxId, 0);
      MeanByFun[ifun].resize(MaxId, 0.);
      StdDevByFun[ifun].resize(MaxId, 0.);
      std::stringstream ii;
      ii << ifun;
      std::string FunFileName = Message::GetResDir() + FullResRootName + ii.str()+ DBext;
      std::ifstream FunDb(FunFileName.c_str(), std::ios_base::in);
      if(!FunDb.is_open())
	{
	  Message::Warning("No database for function %d found... I hope there is no work there because it's gonna be erazed by new results...", ifun);
	  continue;
	}
      else
	{
	  //get the number of points treated for this function
	  int Npoints;
	  std::string line;
	  if(getline (FunDb,line))
	     Npoints = atoi(line.c_str());
	  else
	    Npoints = 0;
	  for (int iP = 0; iP < Npoints; iP++)
	    {
	      //Read some results
	      int id, nres;
	      double x, y, average, stddev;
	      FunDb >> id >> x >> y >> nres >> average >> stddev;
	      NResByFun[ifun][id] = nres;
	      MeanByFun[ifun][id] = average;
	      StdDevByFun[ifun][id] = stddev;
	    }
	  FunDb.close();
	}
    }
}

void Database::ParsePointFiles(){
  Message::Info("ParsePointFiles...");


}

void Database::PrintPointsDone(){
  for (int i = 0; i < NpointsDone; i++)
    PointsDone[i]->Print();
}
