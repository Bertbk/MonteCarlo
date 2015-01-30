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
  Npoints = 0;
  NpointsToDo = 0;
  NResByFun.resize(Message::GetNFUN());
}

Database::~Database(){
  for (int i = 0; i < Npoints; i++)
    delete Points[i];
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
      Npoints = 0;
      return;
    }
  else
    {
      Message::Info("Database file found! Let's rock!");
      //      std::string line;
      //      getline(PointsDb, line);
      PointsDb >> Npoints;
      MaxId = Npoints;
      Points.resize(Npoints);
      for (int i = 0; i < Npoints; i++)
	{
	  int Id;
	  double xi, y;
	  PointsDb >> Id >> xi >> y;
	  Points[i] = new Point(Id, xi, y);
	}
      PointsDb.close();
    }
}

void Database::ParsePointFiles(){
  Message::Info("ParsePointFiles...");
  for (int iPoint = 0; iPoint < Npoints ; iPoint ++)
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
	      //	      FunDbWrite << 0. << std::endl; // Average
	      // FunDbWrite << 0. << std::endl; // Std Dev
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
	      Points[iPoint]->SetNResFiles(ifun, nresfiles);
	      Points[iPoint]->SetMCDone(ifun, nMCDone);
	      //	      Points[iPoint]->SetAverage(ifun, aver);
	      //	      Points[iPoint]->SetStdDev(ifun, standddev);
	    }
	}
    }
}


void Database::UpdatePointsToDo(std::vector<double> Xi, std::vector<double> Y, std::vector<int> MCToDo){
  int nxi = Xi.size();
  int ny = Y.size();
  int nmc = MCToDo.size();

  if(nxi != ny)
    {
      Message::Warning("Vector Xi and Y not of the same size ! Abording...");
      Message::Finalize(EXIT_FAILURE);
    }
  PointsToDo.resize(nxi);
  std::vector<int> newpointindex; //Save indices of points ...
  for (int i = 0; i < nxi; i++)
    {
      double x = Xi[i];
      double y = Y[i];
      //Find database id of the point
      int id = FindPoint(x,y);
      newpointindex.reserve(nxi);
      if(id > -1)
	PointsToDo[i] = Points[i];
      else // build new point
	{
	  id = MaxId;
	  PointsToDo[i] = new Point(id, x, y);
	  newpointindex.push_back(id);
	  MaxId ++;
	}
    }
  Points.resize(MaxId); //Adding new points...
  for (int i =0 ;i < newpointindex.size(); i ++)
    {
      int id = newpointindex[i];
      Points[id] = new Point(id, Xi[i], Y[i]);
      //Build folder
      BuildFolderPoint(id);
      Npoints++;
    }
  //Update points.db
  RebuildPointsDb();
  //Now compare the MC to do/ already done
  
}

void Database::BuildFolderPoint(int id)
{
      std::stringstream iiPoint, backslash;
      backslash  << "/";
      iiPoint << id;
      std::string PointFolder = Message::GetResDir() + PointFolderRootName + iiPoint.str() + backslash.str();
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
	  std::string FunFileName = PointFolder + FunResRootName + iifun.str() + DBext;
	  std::ifstream FunDb(FunFileName.c_str(), std::ios_base::in);
	  if(!FunDb.is_open())
	    {
	      Message::Warning("BuildFolderPoint: %s already exists in %s!!!", FunFileName.c_str(), FunFolder.c_str());
	      FunDb.close();
	    }
	  else{
	    //Create summary file funXX.db
	    std::ofstream FunDbWrite(FunFileName.c_str(), std::ios_base::out);
	    FunDbWrite << 0 << std::endl; // MC done (total)
	    FunDbWrite << 0 << std::endl; // N files
	    FunDbWrite.close();
	  }
	}
}

void Database::RebuildPointsDb()
{
  //Backup file
  std::string backupcommand;
  backupcommand = "cp " + Message::GetResDir() + PointDatabase + DBext + Message::GetResDir() + PointDatabase + DBext + "BACKUP"; 
  system(backupcommand.c_str());
  //Open file: write (eraze)
  std::string PointsDbName = Message::GetResDir() + PointDatabase + DBext;
  std::ofstream Pointsdb(PointsDbName.c_str(), std::ios_base::in);  
  Pointsdb << Npoints << std::endl;
  //Write all points
  for (int i = 0; i < MaxId; i ++)
    {
      Pointsdb << Points[i]->GetId() << Points[i]->GetXi() << Points[i]->GetY() << std::endl;      
    }
  Pointsdb.close();
  //Delete backup file
  backupcommand = "rm " + Message::GetResDir() + PointDatabase + DBext + "BACKUP"; 
  system(backupcommand.c_str());
}


int Database::FindPoint(double x, double y)
{
  return -1;
}

void Database::PrintPoints(){
  for (int i = 0; i < Npoints; i++)
    Points[i]->Print();
}
