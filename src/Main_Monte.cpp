// 2D Monte Carlo Method
// Results are stored in a new directory

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <math.h>
#include <time.h>
#include <vector>

/*#include<mpi.h>
  #include<omp.h>*/

#include "Message.h"
#include "Mesh.h"
#include "Point.h"
#include "Database.h"

using namespace std;

int main(int argc, char *argv[])
{
  //initiatilization (reading arguments, launching MPI,...)
  Message::Initialize(argc, argv);
  //Build database
  Database Db(Message::GetResDir());
  Db.Init();
  Db.PrintPoints();
  //Seed of rand function
  srand(time(NULL) - 360000*Message::GetRank());
  //Reading which points have been done
  if(Message::GetComputeMC())
    {
      if(Message::RootMpi())
	Message::Info("Let's compute some MC...");
      Db.UpdatePointsToDo(Message::GetGridXi(), Message::GetGridY(), Message::GetDesiredMC());
      Db.LaunchMCSimulations();
    }
  if(Message::GetPos())
    {
      Db.PostProcessing(); // write file funXX.pos on root folder
    }
  if(Message::GetGmsh())
    Db.PrintPOS(Message::GetGMSHFileName()); // write .pos files (need GMSH)
  //Exit smoothly like a boss
  Message::Finalize(EXIT_SUCCESS);
  return 0;
}

