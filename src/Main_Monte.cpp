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
      Message::Info("Let's compute some MC...");
      Db.UpdatePointsToDo(Message::GetGridXi(), Message::GetGridY(), Message::GetDesiredMC());
      Db.LaunchMCSimulations();
    }
  
  //compute (or only recompute) average+std deviation only
  //If(Message::GetPos)
  //   Loop on every Done Point
  //   Read file + compute average + standard deviation and store results
  //EndIf

  //If (print)
  //  Build Geo File
  //  Compute Mesh (?) (system(GMSH...)
  //  Re-read mesh to get elements (parser...)
  //  Print Point Done on disk according to GMSH syntaxe and for every elements...
  //EndIf
  
  //Destroy PointDone and PointToDo
  //Exit smoothly
  Message::Finalize(EXIT_SUCCESS);
  return 0;
}

