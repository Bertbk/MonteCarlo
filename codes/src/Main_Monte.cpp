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

using namespace std;

int main(int argc, char *argv[])
{
//initiatilization (reading arguments, launching MPI,...)
Message::Initialize(argc, argv);
//Reading which points have been done
std::vector<Point*> PointDone;
std::vector<Point*> PointToDo;
Point::ReadAllPoints(&PointDone);
if(Message::GetComputeMC) Point::CreatePointsToDo(&PointToDo, &PointDone);

//MPI : decide who does what
// If(Message::GetComputeMC) Then...
//   Loop on every Point
//     Compute the MC simulations
//      Store on disk
//     If(Message::GetPos) Then Also compute Average+Std Deviation and store everything
//     Concatenate files (local to a point)
//   End Loop
//   Concatenate files (summary file)
// EndIf

//Weird stuff but possible: recompute average+std deviation only
//If(Message::GetPos && !Message::GetComputeMC)
//   Loop on every Done Point
//   Read file + compute average + standard deviation and store results
//EndIf

//If (print)
//  Build Geo File
//  Compute Mesh (?) (system(GMSH...)
//  Re-read mesh to get elements (parser...)
//  Print Point on disk according to GMSH syntaxe and for every elements...
//EndIf

//Destroy PointDone and PointToDo
//Exit smoothly
Message::Finalize(EXIT_SUCCESS);
return 0;
}

