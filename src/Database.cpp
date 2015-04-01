#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include <sys/types.h>
#include <sys/stat.h>

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
    {
    Message::Info("Init Database...");
    InitResFolder();
    //Read folder then subfolder, then ...
    int DbExists = ParseRootFiles();
    if(DbExists)
      ParsePointFiles();
    }
  //Broadcase the Database to every other one
  Broadcast(0);
}

void Database::InitResFolder()
{
  if(!Message::RootMpi())
    return;
  //Create result folder (if does not exist)
  mkdir(m_resDir.c_str(), 0700);
  std::string helpDir = m_resDir + Message::GetHelpDir();
  mkdir(helpDir.c_str(), 0700);
  //Check help files
  std::string HelpPointsDbName = Message::GetResDir() + Message::GetHelpDir() +  Message::GetPointDatabase() + Message::GetDBext() + "_help";
  std::ofstream OutHelpPointsDb(HelpPointsDbName.c_str(), std::ios_base::out);
  OutHelpPointsDb << "1: Number of Points" << "\n";
  OutHelpPointsDb << "2: Id_0   Xi_0   Y_0" << "\n";
  OutHelpPointsDb << "3: Id_1   Xi_1   Y_1" << "\n";
  OutHelpPointsDb << "4: Id_2   Xi_2   Y_2" << "\n";
  OutHelpPointsDb << "   .      .      ." << "\n";
  OutHelpPointsDb << "   .      .      ." << "\n";
  OutHelpPointsDb << "   .      .      ." << "\n";
  OutHelpPointsDb.close();

  std::string PointDir = helpDir + Message::GetPointFolderRootName() + "XX/" ;
  mkdir(PointDir.c_str(), 0700);
  std::string HelpFunXXDbName = PointDir + Message::GetFunResRootName() + "XX" + Message::GetDBext() +"_help";
  std::ofstream OutHelpFunXXDb(HelpFunXXDbName.c_str(), std::ios_base::out);
  OutHelpFunXXDb << "1: Number of \""+ Message::GetPointResRootName() + "XX"+Message::GetDBext() +"\" files in funXX folder" << "\n";
  OutHelpFunXXDb << "2: Total number of MC (Monte Carlo simulations)" << "\n";
  /*  OutHelpFunXXDb << "Number of MC in file 0" << "\n";
  OutHelpFunXXDb << "Number of MC in file 1" << "\n";
  OutHelpFunXXDb << "Number of MC in file 2" << "\n";
  OutHelpFunXXDb << "           .          " << "\n";
  OutHelpFunXXDb << "           .          " << "\n";
  OutHelpFunXXDb << "           .          " << "\n";
  OutHelpFunXXDb << "Number of MC in file M-1" << "\n";*/
  OutHelpFunXXDb.close();

  std::string FunDir = PointDir + Message::GetFunResFolderRootName() + "XX/" ;
  mkdir(FunDir.c_str(), 0700);
  std::string HelpPointResXXDbName = FunDir + Message::GetPointResRootName() + "XX" + Message::GetDBext() +"_help";
  std::ofstream OutHelpPointResXXDb(HelpPointResXXDbName.c_str(), std::ios_base::out);
  OutHelpPointResXXDb << "1:   Number of MC (Monte Carlo Simulations) (say M)" << "\n";
  OutHelpPointResXXDb << "2:   Result 0" << "\n";
  OutHelpPointResXXDb << "3:   Result 1" << "\n";
  OutHelpPointResXXDb << "4:   Result 2" << "\n";
  OutHelpPointResXXDb << "        .    " << "\n";
  OutHelpPointResXXDb << "        .    " << "\n";
  OutHelpPointResXXDb << "        .    " << "\n";
  OutHelpPointResXXDb << "M+1: Result M-1" << "\n";
}


int Database::ParseRootFiles(){
  if(!Message::RootMpi())
    return -1;
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
  if(!Message::RootMpi())
    return;
  Message::Info("ParsePointFiles...");
  for (int iPoint = 0; iPoint < Points.size() ; iPoint ++)
    {
      std::stringstream iiPoint;
      iiPoint << iPoint;
      //Create - if not exist - the folder file
      std::string PointFolder = Message::GetResDir() + Message::GetPointFolderRootName() + iiPoint.str() + Message::GetBackSlash();
      mkdir(PointFolder.c_str(), 0700);
	  
      for (int ifun = 0; ifun < Message::GetNFUN() ; ifun ++)
	{
	  std::stringstream iifun;
	  iifun << ifun;
	  //Check if folder exists
	  std::string FunFolder = PointFolder + Message::GetFunResFolderRootName() + iifun.str();
	  mkdir(FunFolder.c_str(), 0700);
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
	      Points[iPoint]->SetNResFiles(ifun, 0);
	      Points[iPoint]->SetMCDone(ifun, 0);
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


void Database::UpdatePointsToDo(std::vector<double> *Xi, std::vector<double> *Y, std::vector<int> *MCToDo){
  if(Message::RootMpi())
    {
      Message::Info("UpdatePointsToDo...");
      int nxi = Xi->size();
      int ny = Y->size();
      int nmc = MCToDo->size();
      Message::Info("nxi = %d ny = %d MaxId = %d Npoints = %d", nxi, ny, Points.size()-1, Points.size());
      if(nxi != ny)
	{
	  Message::Warning("Vector Xi and Y not of the same size ! Abording...");
	  Message::Finalize(EXIT_FAILURE);
	}
      PointsIdToDo.reserve(nxi);
      std::vector<int> newpointindex; //Save indices of points ...
      newpointindex.reserve(nxi);
      for (int i = 0; i < nxi; i++)
	{
	  double xi = (*Xi)[i];
	  double y = (*Y)[i];
	  //Find database id of the point
	  int id = FindPoint(xi,y);
	  if(id > -1)
	    PointsIdToDo.push_back(i); //kn
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
	  //rewrite point database file
	  RebuildPointsDb();
	}
      PreparePointsToDo();
      NpointsToDo = PointsIdToDo.size();
    }
  Broadcast(0);
  //Send id of the points to do...
#ifdef HAVE_MPI
  MPI_Bcast(&NpointsToDo, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if(!Message::RootMpi())
    PointsIdToDo.resize(NpointsToDo);
  MPI_Bcast(&PointsIdToDo[0], NpointsToDo, MPI_INT, 0, MPI_COMM_WORLD);
#endif
}

void Database::BuildFolderPoint(int id)
{
  if(!Message::RootMpi())
    return;  
  std::stringstream iiPoint;
  iiPoint << id;
  std::string PointFolder = Message::GetResDir() + Message::GetPointFolderRootName() + iiPoint.str() + Message::GetBackSlash();
  //Create - if not exist - the folder file
  mkdir(PointFolder.c_str(), 0700);
  for (int ifun = 0; ifun < Message::GetNFUN() ; ifun ++)
    {
      std::stringstream iifun;
      iifun << ifun;
      //Check if folder exists
      std::string FunFolder = PointFolder + Message::GetFunResFolderRootName() + iifun.str();
      mkdir(FunFolder.c_str(), 0700);
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
	FunDbWrite << 0 << std::endl; // N files
	FunDbWrite << 0 << std::endl; // MC done (total)
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
  //Open file: write (eraze)
  std::string PointsDbName = Message::GetResDir() + Message::GetPointDatabase() + Message::GetDBext();
  std::ofstream Pointsdb(PointsDbName.c_str(), std::ios_base::out);
  Pointsdb << Points.size() << std::endl;
  //Write all points
  for (int i = 0; i < Points.size(); i ++)
      Pointsdb << Points[i]->GetId() << " " << Points[i]->GetXi() << " " << Points[i]->GetY() << std::endl;      
  Pointsdb.close();
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
  double precision = 1e-8;
  for (int i = 0; i < Points.size(); i++)
    {
      double thisXi = Points[i]->GetXi();
      double thisY = Points[i]->GetY();
      if(abs(thisXi- xi)<=precision && abs(thisY-y)<=precision)
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
  int nPToPos = Points.size();
  if(nPToPos == 0)
    {
      Message::Warning("Post processing impossible: no point detected");
      return;
    }
  for (int ifun =0; ifun < NFUN; ifun ++)
    {
      //Folder name
      std::stringstream iifun;
      iifun << ifun;
      std::string rootFunFolder = Message::GetFunResFolderRootName() + iifun.str() + Message::GetBackSlash();
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

void Database::PrintGMSHPOS()
{
  if(!Message::RootMpi())
    return;
  int np = Points.size();
  if(np == 0)
    {
      Message::Warning("GMSH Post processing impossible: no point detected");
      return;
    }
  //Rebuild Database to be sure to have the last informations
  std::vector<double> Xi;
  std::vector<double> Y;
  std::vector<std::vector<double> > res;

  ExtractXYRes(&Xi, &Y, &res);
  Mesh myMesh(Xi, Y);
  myMesh.SetRes(&res);
  myMesh.Update();
  myMesh.PrintGMSHRes(Message::GetResDir() + Message::GetFunResRootName());
  myMesh.PrintMesh(Message::GetResDir()+"mesh");
}

void Database::Broadcast(int rank)
{
#ifdef HAVE_MPI
  if(Message::GetRank() == rank)
    Message::Info("Sending database");
  else
    Message::Info("Receiving database from rank %d", rank);
  int Npoints;
  if(Message::GetRank()==rank)
    Npoints = Points.size();
  //Bcast Npoint
  MPI_Bcast(&Npoints, 1, MPI_INT, rank, MPI_COMM_WORLD);
  if(Message::GetRank() !=rank)
    {
      //Clean the vectors of eventual preceeding points
      for (int iP = 0; iP < Points.size(); iP++)
	delete Points[iP];
      Points.resize(Npoints);
      for (int iP = 0; iP < Npoints; iP++)
	Points[iP] = NULL;
    }
  //Broadcast every point
  for (int iP = 0; iP < Npoints; iP++)
    {
      int id, xi, y;
      if(Message::GetRank() ==rank)
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
