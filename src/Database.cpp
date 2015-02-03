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
  NpointsToDo = 0;
  //  NResByFun.resize(Message::GetNFUN());
}

Database::~Database(){
  for (int i = 0; i < Points.size(); i++)
    delete Points[i];
  for (int i = 0; i < PointsToDo.size(); i++)
    delete PointsToDo[i];
}

void Database::Init()
{
  Message::Info("Init Database...");
  //Read folder then subfolder, then ...
  int DbExists = ParseRootFiles();
  if(DbExists)
    ParsePointFiles();

}


int Database::ParseRootFiles(){
  Message::Info("ParseRootFiles...");
  //The root directory should have a file points.db and N_FUN files res_funXX.db
  //points.db: contains information about the points (Id, X, Y)
  std::string PointsDbName = Message::GetResDir() + PointDatabase + DBext;
  std::ifstream PointsDb(PointsDbName.c_str(), std::ios_base::in);
  if(!PointsDb.is_open())
    {
      Message::Warning("No database file found... I hope there is no work there because it's gonna be erazed by new results...");
      return 0;
    }
  else
    {
      Message::Info("Database file found! Let's rock!");
      //      std::string line;
      //      getline(PointsDb, line);
      int Npoints;
      PointsDb >> Npoints;
      Points.resize(Npoints);
      for (int i = 0; i < Npoints; i++)
	{
	  int Id;
	  double xi, y;
	  PointsDb >> Id >> xi >> y;
	  Points[i] = new Point(Id, xi, y);
	}
      PointsDb.close();
      return 1;
    }
}

void Database::ParsePointFiles(){
  Message::Info("ParsePointFiles...");
  for (int iPoint = 0; iPoint < Points.size() ; iPoint ++)
    {
      std::stringstream iiPoint, backslash;
      backslash  << "/";
      iiPoint << iPoint;
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
	  //Read summary file funXX.db
	  std::string FunFileName = PointFolder + FunResRootName + iifun.str() + DBext;
	  //Open file
	  std::ifstream FunDb(FunFileName.c_str(), std::ios_base::in);
	  if(!FunDb.is_open())
	    {
	      Message::Warning("No %s for Point %d and function %d found... Building empty file", FunResRootName.c_str(), iPoint, ifun);
	      std::ofstream FunDbWrite(FunFileName.c_str(), std::ios_base::out);
	      FunDbWrite << 0 << std::endl; // MC done (total)
	      FunDbWrite << 0 << std::endl; // N files
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
	    }
	}
    }
  Message::Info("End ParsePointFiles...");
}


void Database::UpdatePointsToDo(std::vector<double> *Xi, std::vector<double> *Y, std::vector<int> *MCToDo){
  Message::Info("UpdatePointsToDo...");
  int nxi = Xi->size();
  int ny = Y->size();
  int nmc = MCToDo->size();
  Message::Info("nxi = %d ny = %d MaxId = %d Npoints = %d", nxi, ny, Points.size(), Points.size());

  if(nxi != ny)
    {
      Message::Warning("Vector Xi and Y not of the same size ! Abording...");
      Message::Finalize(EXIT_FAILURE);
    }
  PointsIdToDo.reserve(nxi);
  std::vector<int> newpointindex; //Save indices of points ...
  for (int i = 0; i < nxi; i++)
    {
      double xi = (*Xi)[i];
      double y = (*Y)[i];
      //Find database id of the point
      int id = FindPoint(xi,y);
      newpointindex.reserve(nxi);
      if(id > -1)
	{
	PointsIdToDo.push_back(i);
	Message::Info("Adding Id %d to PointsIdToDo", i);
	}
      else // build new point
	{
	  //They will be created and added below, after the loop
	  newpointindex.push_back(i);
	}
    }
  int newId = Points.size();
  Points.resize(Points.size() + newpointindex.size()); //Adding new points...
  if(newpointindex.size()>0)
    {
      Message::Info("Creating new folders and rebuilding the Database Points.db");
      for (int i =0 ;i < newpointindex.size(); i ++)
	{
	  int ixi = newpointindex[i];
	  Points[newId] = new Point(newId, (*Xi)[ixi], (*Y)[ixi]);
	  PointsIdToDo.push_back(newId);
	  newId++;
	  //Build folder
	  BuildFolderPoint(newId);
	}
      //Update points.db
      RebuildPointsDb();
    }
  PreparePointsToDo();
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
	  if(FunDb.is_open())
	    {
	      Message::Warning("BuildFolderPoint: %s already exists in %s!!!", FunFileName.c_str(), PointFolder.c_str());
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
  Message::Info("RebuildPointsDb...");
  //Backup file
  std::string backupcommand;
  backupcommand = "cp " + Message::GetResDir() + PointDatabase + DBext + Message::GetResDir() + PointDatabase + DBext + "BACKUP"; 
  system(backupcommand.c_str());
  //Open file: write (eraze)
  std::string PointsDbName = Message::GetResDir() + PointDatabase + DBext;
  std::ofstream Pointsdb(PointsDbName.c_str(), std::ios_base::out);
  Pointsdb << Points.size() << std::endl;
  //Write all points
  for (int i = 0; i < Points.size(); i ++)
      Pointsdb << Points[i]->GetId() << " " << Points[i]->GetXi() << " " << Points[i]->GetY() << std::endl;      
  Pointsdb.close();
  //Delete backup file
  backupcommand = "rm " + Message::GetResDir() + PointDatabase + DBext + "BACKUP"; 
  system(backupcommand.c_str());
}

void Database::PreparePointsToDo()
{
  Message::Info("PreparePointsToDo...");
  int npToDo = PointsIdToDo.size();
  std::vector<int> *desiredMC = Message::GetDesiredMC();
  for (int iP =0 ; iP < npToDo; iP++)
    {
      int id = PointsIdToDo[iP];
      Point *myP = Points[id];
      myP->SetMCToDo(desiredMC);
    }
}



int Database::FindPoint(double xi, double y)
{
  for (int i = 0; i < Points.size(); i++)
    {
      double thisXi = Points[i]->GetXi();
      double thisY = Points[i]->GetY();
      if(thisXi == xi && thisY == y)
	return Points[i]->GetId();
    }
  return -1;
}

void Database::PrintPoints(){
  Message::Info("PrintPoints...");
  for (int i = 0; i < Points.size(); i++)
    Points[i]->Print();
}
