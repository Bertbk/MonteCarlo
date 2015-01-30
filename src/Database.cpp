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
std::string Database::PointFolderRootName = "Point";
std::string Database::FunResFolderRootName = "fun";
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
  ParsePointFiles();

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
  for (int iPoint = 0; iPoint < NpointsDone ; iPoint ++)
    {
      std::stringstream iiPoint, backslash;
      backslash  << "/";
      iiPoint << iPoint;
      std::string PointFolder = Message::GetResDir() + PointFolderRootName + iiPoint.str() + backslash.str();
      Message::Debug("PointFolder : %s", PointFolder.c_str());
      //Create - if not exist - the foler file
      std::string command_PointFolder = "if ! test -d " + PointFolder + "; then mkdir "+ PointFolder+"; fi";
      system(command_PointFolder.c_str());

      for (int ifun = 0; ifun < Message::GetNFUN() ; ifun ++)
	{
	  std::stringstream iifun;
	  iifun << ifun;
	  //Check if folder exists
	  std::string FunFolder = PointFolder + FunResFolderRootName + iifun.str();
	  std::string command_FunFolder = "if ! test -d " + FunFolder + "; then mkdir "+ FunFolder+"; fi";
	  system(command_FunFolder.c_str());
	  //Read summary file funXX.db
	  std::string FunFileName = PointFolder + FunResRootName + iifun.str() + DBext;
	  Message::Debug("FunFileName : %s", FunFileName.c_str());
	  //Open file
	  std::ifstream FunDb(FunFileName.c_str(), std::ios_base::in);
	  if(!FunDb.is_open())
	    {
	      Message::Warning("No %s for Point %d and function %d found... Building empty file", FunResRootName.c_str(), iPoint, ifun);
	      std::ofstream FunDbWrite(FunFileName.c_str(), std::ios_base::out);
	      FunDbWrite << 0 << std::endl; // MC done (total)
	      FunDbWrite << 0 << std::endl; // N files
	      FunDbWrite << 0. << std::endl; // Average
	      FunDbWrite << 0. << std::endl; // Std Dev
	      FunDbWrite.close();
	      FunDb.open(FunFileName.c_str(), std::ios_base::in);
	    }
	  else
	    {
	      //Read files
	      int nresfiles, nMCDone;
	      double aver, standddev;
	      FunDb >> nresfiles;
	      FunDb >> nMCDone;
	      FunDb >> aver;
	      FunDb >> standddev;
	      PointsDone[iPoint]->SetNResFiles(ifun, nresfiles);
	      PointsDone[iPoint]->SetMCDone(ifun, nMCDone);
	      PointsDone[iPoint]->SetAverage(ifun, aver);
	      PointsDone[iPoint]->SetStdDev(ifun, standddev);
	    }
	}
    } 
}

void Database::PrintPointsDone(){
  for (int i = 0; i < NpointsDone; i++)
    PointsDone[i]->Print();
}
