#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "Database.h"
#include "Message.h"
#include "Point.h"

std::string Database::DBext = Message::GetDBext();
std::string Database::POSext = Message::GetPOSext();
std::string Database::PointDatabase = Message::GetPointDatabase();
std::string Database::FullResRootName = Message::GetFullResRootName();
std::string Database::CurrentPointDatabase = Message::GetCurrentPointDatabase();
std::string Database::PointFolderRootName = Message::GetPointFolderRootName();
std::string Database::FunResFolderRootName = Message::GetFunResFolderRootName();
std::string Database::FunResRootName = Message::GetFunResRootName();
std::string Database::PointResRootName = Message::GetPointResRootName();
std::string Database::BackSlash = Message::GetBackSlash();

//Constructor
Database::Database(std::string resdir){
  m_resDir = resdir;
  NpointsToDo = 0;
}

Database::~Database(){
  for (int i = 0; i < Points.size(); i++)
    delete Points[i];
}

void Database::Init()
{
  Message::Info("Init Database...");
  CheckOrBuildRootFolder();
  //Read folder then subfolder, then ...
  int DbExists = ParseRootFiles();
  if(DbExists)
    {
      ParsePointFiles();
      ParsePosFiles();
    }
}

void Database::CheckOrBuildRootFolder()
{
  if(!Message::RootMpi())
    return;
  //Create result folder (if does not exist)
  std::string command = "if ! test -d " + m_resDir + "; then mkdir "+ m_resDir+"; fi";
  system(command.c_str());
  std::string helpDir = m_resDir + Message::GetHelpDir();
  command = "if ! test -d " + helpDir + "; then mkdir "+ helpDir+"; fi";
  system(command.c_str());
  //Check help files
  std::string HelpPointsDbName = Message::GetResDir() + Message::GetHelpDir() +  Message::GetPointDatabase() + Message::GetDBext() + "_help";
  std::ofstream OutHelpPointsDb(HelpPointsDbName.c_str(), std::ios_base::out);
  OutHelpPointsDb << "Number of Points" << "\n";
  OutHelpPointsDb << "Id_0   Xi_0   Y_0" << "\n";
  OutHelpPointsDb << "Id_1   Xi_1   Y_1" << "\n";
  OutHelpPointsDb << "Id_2   Xi_2   Y_2" << "\n";
  OutHelpPointsDb << " .      .      ." << "\n";
  OutHelpPointsDb << " .      .      ." << "\n";
  OutHelpPointsDb << " .      .      ." << "\n";
  OutHelpPointsDb.close();

  std::string PointDir = helpDir + Message::GetPointFolderRootName() + "XX/" ;
  command = "if ! test -d " + PointDir + "; then mkdir "+ PointDir +"; fi";
  system(command.c_str());
  std::string HelpFunXXDbName = PointDir + Message::GetFunResRootName() + "XX" + Message::GetDBext() +"_help";
  std::ofstream OutHelpFunXXDb(HelpFunXXDbName.c_str(), std::ios_base::out);
  OutHelpFunXXDb << "Number of \""+ Message::GetPointResRootName() + "XX"+Message::GetDBext() +"\" files in funXX folder (say M)" << "\n";
  OutHelpFunXXDb << "Total number of MC (Monte Carlo simulations)" << "\n";
  OutHelpFunXXDb << "Number of MC in file 0" << "\n";
  OutHelpFunXXDb << "Number of MC in file 1" << "\n";
  OutHelpFunXXDb << "Number of MC in file 2" << "\n";
  OutHelpFunXXDb << "           .          " << "\n";
  OutHelpFunXXDb << "           .          " << "\n";
  OutHelpFunXXDb << "           .          " << "\n";
  OutHelpFunXXDb << "Number of MC in file M-1" << "\n";
  OutHelpFunXXDb.close();

  std::string FunDir = PointDir + Message::GetFunResFolderRootName() + "XX/" ;
  command = "if ! test -d " + FunDir + "; then mkdir "+ FunDir +"; fi";
  system(command.c_str());
  std::string HelpPointResXXDbName = FunDir + Message::GetPointResRootName() + "XX" + Message::GetDBext() +"_help";
  std::ofstream OutHelpPointResXXDb(HelpPointResXXDbName.c_str(), std::ios_base::out);
  OutHelpPointResXXDb << "Number of MC (Monte Carlo Simulations) (say M)" << "\n";
  OutHelpPointResXXDb << "Result 0" << "\n";
  OutHelpPointResXXDb << "Result 1" << "\n";
  OutHelpPointResXXDb << "Result 2" << "\n";
  OutHelpPointResXXDb << "   .    " << "\n";
  OutHelpPointResXXDb << "   .    " << "\n";
  OutHelpPointResXXDb << "   .    " << "\n";
  OutHelpPointResXXDb << "Result M-1" << "\n";
}


int Database::ParseRootFiles(){
  Message::Info("ParseRootFiles...");
  //The root directory should have a file points.db and N_FUN files res_funXX.db
  //points.db: contains information about the points (Id, X, Y)
  std::string PointsDbName = Message::GetResDir() + PointDatabase + DBext;
  std::ifstream PointsDb(PointsDbName.c_str(), std::ios_base::in);
  if(!PointsDb.is_open())
    {
      if(Message::RootMpi())
	Message::Warning("No database file found... I hope there is no work there because it's gonna be erazed by new results...");
      return 0;
    }
  else
    {
      Message::Info("Database file found! Let's rock!");
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

void Database::ParsePosFiles(){
  Message::Info("ParsePosFiles...");
}

void Database::ParsePointFiles(){
  Message::Info("ParsePointFiles...");
  for (int iPoint = 0; iPoint < Points.size() ; iPoint ++)
    {
      std::stringstream iiPoint;
      iiPoint << iPoint;
      std::string PointFolder = Message::GetResDir() + PointFolderRootName + iiPoint.str() + BackSlash;
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
	  //Build folder
	  BuildFolderPoint(newId);
	  //Increase NewId size
	  newId++;
	}

      //Update points.db
#if defined HAVE_MPI
      //Possible problem with reading process still searching the file...
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      RebuildPointsDb();
    }
  PreparePointsToDo();
}

void Database::BuildFolderPoint(int id)
{
  if(!Message::RootMpi())
    return;  
  std::stringstream iiPoint;
  iiPoint << id;
  std::string PointFolder = Message::GetResDir() + PointFolderRootName + iiPoint.str() + BackSlash;
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
  if(!Message::RootMpi())
    return;
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

void Database::LaunchMCSimulations()
{
  Message::Info("LaunchMCSimulations...");
  int npToDo = PointsIdToDo.size();
  int myRank = Message::GetRank();
  std::vector<int> iP_start, iP_end;
  Message::DistributeWork(npToDo, &iP_start, &iP_end);
  for (int iP = iP_start[myRank] ; iP < iP_end[myRank]; iP++)
    {
      int id = PointsIdToDo[iP];
      Points[id]->LaunchMC();
    }
#ifdef HAVE_MPI
  //Exchange information about the point
  int nMpi = Message::GetNProc();
  for (int iRank = 0 ; iRank < nMpi; iRank++)
    {
      for (int iP = iP_start[iRank] ; iP < iP_end[iRank]; iP++)
	{
	  int id = PointsIdToDo[iP];
	  Points[id]->Broadcast(iRank);
	}
    }
#endif
}


void Database::PrintPoints(){
  if(!Message::RootMpi())
    return;
  Message::Info("PrintPoints...");
  for (int i = 0; i < Points.size(); i++)
    Points[i]->Print();
}

void Database::PostProcessing(){
  int NFUN = Message::GetNFUN();
  for (int ifun =0; ifun < NFUN; ifun ++)
    {
      //Folder name
      std::stringstream iifun;
      iifun << ifun;
      std::string rootFunFolder = FunResFolderRootName + iifun.str() + BackSlash;
      //This loop is parallelizable (MPI)
      int nPToPos = Points.size();
      int myRank = Message::GetRank();
      int nMpi = Message::GetNProc();
      std::vector<int> iP_start, iP_end;
      Message::DistributeWork(nPToPos, &iP_start, &iP_end);
      for (int iP = iP_start[myRank] ; iP < iP_end[myRank]; iP++)
	  Points[iP]->PostProcessing(ifun);
#if defined HAVE_MPI
      //Exchange of information (give everything to proc 0)
      std::vector<MPI_Request> request(0);
      if(myRank == 0)
	{
	  for (int iRank = 1; iRank < nMpi; iRank ++)
	    {
	      for (int iP = iP_start[iRank] ; iP < iP_end[iRank]; iP++)
		Points[iP]->Irecv(iRank, &request);
	    }
	}
      else
	{
	  for (int iP = iP_start[myRank] ; iP < iP_end[myRank]; iP++)
	    Points[iP]->Isend(0, &request);
	}
      std::vector< MPI_Status > tab_status(request.size());
      MPI_Waitall(request.size(), &request[0], &tab_status[0]);
#endif
      if(Message::RootMpi())
	{
	  //Write on files !
	  std::string funXXName = Message::GetResDir() + FunResRootName +iifun.str() + POSext;
	  std::ofstream funXX(funXXName.c_str(), std::ios_base::out);
	  funXX << Points.size() << "\n";
	  for (int iP = 0; iP < Points.size(); iP++)
	    funXX << Points[iP]->GetXi() << " "<< Points[iP]->GetY() << " "<< Points[iP]->GetAverage(ifun) << " "<< Points[iP]->GetStdDev(ifun) << " " << "\n";
	  funXX.close();
	}
    }
}


void Database::PostProcessingGMSH()
{
  //Rebuild Database to be sure to have the last informations
  
  
  
  
}
