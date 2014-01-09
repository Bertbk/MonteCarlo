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
//#include "Point.h"

using namespace std;

int main(int argc, char *argv[])
{
  //initiatilization (reading arguments, launching MPI,...)
  Message::Initialize(argc, argv);
  Message::Finalize(EXIT_SUCCESS);
  return 0;
}

