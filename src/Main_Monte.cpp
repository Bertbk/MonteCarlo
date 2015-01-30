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
  //MPI : decide who does what ...
  std::vector<Point*> PointDone;
  Point::ReadAllPoints(&PointDone);
  
  if(Message::GetComputeMC())
    {
      std::vector<Point*> PointToDo;
      Point::CreatePointsToDo(&PointToDo, &PointDone);
      int nPointToDo=PointToDo.size();
      std::vector<int> IndexOfPointToDo; // in sequential, this will be 0:(nPointToDo-1)
      Message::DistributeWork(nPointToDo, &IndexOfPointToDo);
      for(int i =0; i < nPointToDo; i++)
	{
	  Point *cPoint = PointToDo[IndexOfPointToDo[i]];
	  //Prepare folder, files,...
	  cPoint->LaunchMC();
	  //    For each MC simulations:
	  //    -- Compute res
	  //    -- Store on disk (resDir/idXX/res_aux file)
	  //    //UNSURE IF I DON'T SPLIT EVERYTHING (SIMPLER)If(Message::GetPos) Then Also compute Average+Std Deviation and store everything
	  //    Concatenate resDir/idXX/res_aux files
	  //    Update DBB file for this point (resDir/dbb_aux)
	}//  End Loop
      // Backup resDir/ddb to resDir/ddb_backup($TIME)
      //  Concatenate resDir/ddb_aux to resDir/ddb (whatever the order of id...)
      //Delete Every points PointDone, and do PointDone = PointToDo
    }//EndIf
  
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

