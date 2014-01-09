// 2D Monte Carlo Method
// Results are stored in a new directory

#include<iostream>
#include<cstdlib>
#include<ctime>
#include<cmath>
#include<fstream>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>
#include <vector>
#include <iomanip>
#include <sstream>


/*#include<mpi.h>
  #include<omp.h>*/

#include "Message.h"
#include "MyResults.h"
#include "Point.h"

using namespace std;

int main(int argc, char *argv[])
{
  //initiatilization (reading arguments, launching MPI,...)
  Message::Initialize(argc, argv);
  Message::Finalize();
  return 0;
}

