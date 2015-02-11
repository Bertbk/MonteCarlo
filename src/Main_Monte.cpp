// 2D Monte Carlo Method
// Results are stored in a directory
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>

#include "MonteCarlo/Message.h"
#include "MonteCarlo/Mesh.h"
#include "MonteCarlo/Point.h"
#include "MonteCarlo/Database.h"


int main(int argc, char *argv[])
{
  //initiatilization (reading arguments, launching MPI,...)
  Message::Initialize(argc, argv);
  //Build database
  Database Db(Message::GetResDir());
  Db.Init();
  Db.PrintPoints();
  if(Message::GetComputeMC())
    {
      //The computation of the monte carlo simulations will be done here
      if(Message::RootMpi())
	Message::Info("Let's compute some MC...");
      Db.UpdatePointsToDo(Message::GetGridXi(), Message::GetGridY(), Message::GetDesiredMC());
      Db.LaunchMCSimulations();
    }
  if(Message::GetPos())
    {
      if(Message::RootMpi())
	Message::Info("Let's compute some MC...");
      // write file funXX.pos on root folder
      Db.PostProcessing(); 
    }
  if(Message::GetGmsh())
    {
      // write .pos files (need GMSH)
      Db.PrintPOS(Message::GetGMSHFileName());
    }
  //Exit smoothly like a boss
  Message::Finalize(EXIT_SUCCESS);
  return 0;
}

