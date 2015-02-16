#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include <MonteCarlo/Database.h>
#include <MonteCarlo/Message.h>
#include <MonteCarlo/Mesh.h>
#include <MonteCarlo/Point.h>

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
  if(Message::RootMpi())
    Message::Info("Init Database...");
  InitResFolder();
  Message::Barrier();
  //Read folder then subfolder, then ...
  int DbExists = ParseRootFiles();
  Message::Barrier();
  if(DbExists)
      ParsePointFiles();
}

void Database::InitResFolder()
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
  std::string PointsDbName = Message::GetResDir() + Message::GetPointDatabase() + Message::GetDBext();
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


void Database::ParsePointFiles(){
  if(Message::RootMpi())
    {
      Message::Info("ParsePointFiles...");
      for (int iPoint = 0; iPoint < Points.size() ; iPoint ++)
	{
	  std::stringstream iiPoint;
	  iiPoint << iPoint;
	  //Create - if not exist - the folder file
	  std::string PointFolder = Message::GetResDir() + Message::GetPointFolderRootName() + iiPoint.str() + Message::GetBackSlash();
	  std::string command_PointFolder = "if ! test -d " + PointFolder + "; then mkdir "+ PointFolder+"; fi";
	  system(command_PointFolder.c_str());
	  
	  for (int ifun = 0; ifun < Message::GetNFUN() ; ifun ++)
	    {
	      std::stringstream iifun;
	      iifun << ifun;
	      //Check if folder exists
	      std::string FunFolder = PointFolder + Message::GetFunResFolderRootName() + iifun.str();
	      std::string command_FunFolder = "if ! test -d " + FunFolder + "; then mkdir "+ FunFolder+"; fi";
	      system(command_FunFolder.c_str());
	      //Read summary file funXX.db
	      std::string FunFileName = PointFolder + Message::GetFunResRootName() + iifun.str() + Message::GetDBext();
	      //Open file
	      std::ifstream FunDb(FunFileName.c_str(), std::ios_base::in);
	      if(!FunDb.is_open())
		{
		  Message::Warning("No %s for Point %d and function %d found... Building empty file", Message::GetFunResRootName().c_str(), iPoint, ifun);
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
		  FunDb >> nresfiles;
		  FunDb >> nMCDone;
		  Points[iPoint]->SetNResFiles(ifun, nresfiles);
		  Points[iPoint]->SetMCDone(ifun, nMCDone);
		}
	    }
	}
      Message::Info("End ParsePointFiles...");
    }
  //Broadcase the Database to every other one
  Broadcast(0);

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
      Message::Info("Building new points (folder, etc.)");
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
      RebuildPointsDb();
    }
  PreparePointsToDo();
  Message::Barrier();
}

void Database::BuildFolderPoint(int id)
{
  if(!Message::RootMpi())
    return;  
  std::stringstream iiPoint;
  iiPoint << id;
  std::string PointFolder = Message::GetResDir() + Message::GetPointFolderRootName() + iiPoint.str() + Message::GetBackSlash();
  //Create - if not exist - the foler file
  std::string command_PointFolder = "if ! test -d " + PointFolder + "; then mkdir "+ PointFolder+"; fi";
  system(command_PointFolder.c_str());
  for (int ifun = 0; ifun < Message::GetNFUN() ; ifun ++)
    {
      std::stringstream iifun;
      iifun << ifun;
      //Check if folder exists
      std::string FunFolder = PointFolder + Message::GetFunResFolderRootName() + iifun.str();
      std::string command_FunFolder = "if ! test -d " + FunFolder + "; then mkdir "+ FunFolder+"; fi";
      system(command_FunFolder.c_str());
      std::string FunFileName = PointFolder + Message::GetFunResRootName() + iifun.str() + Message::GetDBext();
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
  backupcommand = "cp " + Message::GetResDir() + Message::GetPointDatabase() + Message::GetDBext() + Message::GetResDir() + Message::GetPointDatabase() + Message::GetDBext() + "BACKUP"; 
  system(backupcommand.c_str());
  //Open file: write (eraze)
  std::string PointsDbName = Message::GetResDir() + Message::GetPointDatabase() + Message::GetDBext();
  std::ofstream Pointsdb(PointsDbName.c_str(), std::ios_base::out);
  Pointsdb << Points.size() << std::endl;
  //Write all points
  for (int i = 0; i < Points.size(); i ++)
      Pointsdb << Points[i]->GetId() << " " << Points[i]->GetXi() << " " << Points[i]->GetY() << std::endl;      
  Pointsdb.close();
  //Delete backup file
  backupcommand = "rm " + Message::GetResDir() + Message::GetPointDatabase() + Message::GetDBext() + "BACKUP"; 
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
  if(Message::RootMpi())
    Message::Info("Post-Processing...");
  int NFUN = Message::GetNFUN();
  for (int ifun =0; ifun < NFUN; ifun ++)
    {
      //Folder name
      std::stringstream iifun;
      iifun << ifun;
      std::string rootFunFolder = Message::GetFunResFolderRootName() + iifun.str() + Message::GetBackSlash();
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
	  std::string funXXName = Message::GetResDir() + Message::GetFunResRootName() +iifun.str() + Message::GetPOSext();
	  std::ofstream funXX(funXXName.c_str(), std::ios_base::out);
	  funXX << Points.size() << "\n";
	  for (int iP = 0; iP < Points.size(); iP++)
	    funXX << Points[iP]->GetXi() << " "<< Points[iP]->GetY() << " "<< Points[iP]->GetAverage(ifun) << " "<< Points[iP]->GetStdDev(ifun) << " " << "\n";
	  funXX.close();
	}
    }
}


void Database::ExtractXYRes(std::vector<double> *Xi,std::vector<double> *Y,std::vector<std::vector<double> > *res)
{
  int np = Points.size();
  Xi->resize(np);
  Y->resize(np);
  res->resize(Message::GetNFUN());
  for (int ifun = 0; ifun < Message::GetNFUN(); ifun++)
    (*res)[ifun].resize(np);
  for (int iP=0; iP < np; iP++)
    {
      Point *p = Points[iP];
      (*Xi)[iP] = p->GetXi();
      (*Y)[iP] = p->GetY();
      for (int ifun = 0; ifun < Message::GetNFUN(); ifun++)
	(*res)[ifun][iP] = p->GetAverage(ifun);
    }
}

void Database::PrintPOS(std::string FileName)
{
  if(!Message::RootMpi())
    return;
  //Rebuild Database to be sure to have the last informations
  std::vector<double> Xi;
  std::vector<double> Y;
  std::vector<std::vector<double> > res;

  ExtractXYRes(&Xi, &Y, &res);
  Mesh myMesh(Xi, Y);
  myMesh.SetRes(&res);
  myMesh.Update();
  myMesh.PrintRes(Message::GetResDir()+FileName);
  myMesh.PrintMesh(Message::GetResDir()+"lol");
}

void Database::Broadcast(int rank)
{
#ifdef HAVE_MPI
  int Npoints;
  if(Message::GetRank()==rank)
    Npoints = Points.size();
  //Bcast Npoint
  MPI_Bcast(&Npoints, 1, MPI_INT, rank, MPI_COMM_WORLD);
  if(Message::GetRank() !=rank)
    Points.resize(Npoints);
  //Broadcast every point
  for (int iP = 0; iP < Npoints; iP++)
    {
      int id, xi, y;
      if(Message::GetRank() !=rank)
	{
	  id = Points[iP]->GetId();
	  xi = Points[iP]->GetXi();
	  y = Points[iP]->GetY();
	}
      MPI_Bcast(&id, 1, MPI_INT, rank, MPI_COMM_WORLD);
      MPI_Bcast(&xi, 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);
      MPI_Bcast(&y, 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);
      if(Message::GetRank() !=rank)
	Points[iP] = new Point(id, xi, y);
      Points[iP]->Broadcast(rank);
    }
#endif
}
