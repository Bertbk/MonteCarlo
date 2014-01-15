#include <iostream>
#include <fstream>
#include <sstream> //for osstream
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <time.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <string>
#include <algorithm> //for remove_if, isspace

#include "Message.h"

#if defined(WITH_MPI)
#include "mpi.h"
#endif

#if defined(WITH_OMP)
#include "omp.h"
#endif


//using namespace std;
int Message::_ComputeMC = 0;
int Message::_Pos = 0;
int Message::_Gmsh = 0;
int Message::_verbosity = 4;
int Message::_myRank = 0;
int Message::_nb_proc = 1;
//PARAMETERS
//==========
double Message::_deuxpi = 8*atan(1.0);
double Message::_lambda = 0.1;
//parameters of the oscillator
double Message::_Y   = 1.0;
double Message::_c0  = 1.0;
double Message::_k   = 1.0;
//filtre
double Message::_alpha = 0.5;  
//correlation of the noises between y and gamma (g)
double Message::_ro = 0.1;
double Message::_roc = sqrt(1.0 - _ro*_ro);
//parameters of the numerical procedure
double Message::_T   = 500.0;
double Message::_dt  = 0.0001;
double Message::_sdt = sqrt(_dt);
//Grid
double Message::_xi_min=-5, Message::_xi_max=5, Message::_dxi = 0.1;
double Message::_dy=0.1, Message::_y_min=0, Message::_y_max=5;
//choice of function (for final computation)
const int Message::_NFUN = 4;
std::vector<int> Message::_FunChoice(_NFUN);
std::vector<int> Message::_desired_MC(_NFUN);
std::string Message::_paramFile = "param";
std::string Message::_resDir = "res/";

//Message
//-------
void Message::Initialize(int argc, char *argv[])
{
#if defined(WITH_MPI)
  MPI_Init( &argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &_myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &_nb_proc);
  Message::Info("Launched with MPI (%d processes)", _nb_proc);
#endif
#if defined(WITH_OMP)
#pragma omp parallel
  {
  _nb_threads = omp_get_num_threads();
  }
  Message::Info("Launched with OpenMP (%d threads)", _nb_threads);
#endif
  int i = 1;
  bool doCheckOnly = 0;
  while (i < argc) {
    if (argv[i][0] == '-') {
      if (!strcmp(argv[i]+1, "check"))    { doCheckOnly = 1; i++;}
      else if (!strcmp(argv[i]+1, "v"))   { _verbosity = atof(argv[i+1]) ; i+=2 ; }
      else if (!strcmp(argv[i]+1, "par")) { _paramFile = argv[i+1]; i+=2; }
      else if (!strcmp(argv[i]+1, "MC"))  { _ComputeMC = 1; i++; }
      else if (!strcmp(argv[i]+1, "pos")) { _Pos = 1; i++; }
    }
    else
      {
	printf("What the hell is this option ?\n");
	i++;
      }
  }
  for(int i =0; i<_NFUN; i++)
    {
      _desired_MC[i] = 0;
      _FunChoice[i] = 0;
    }
  //Parse param file
  Message::Parse();
  //Print info
  if(doCheckOnly) { Message::Check(); Message::Finalize(EXIT_SUCCESS);}
  else {if (_verbosity > 0) Message::Check();}

  //Create result folder (if does not exist)
  std::string command = "if ! test -d " + _resDir + "; then mkdir "+ _resDir+"; fi";
  system(command.c_str());
}


//Info...
void Message::Info(const char *format, ...)
{
  char str[1024];
  va_list args;
  va_start (args, format);
  vsnprintf (str, 1024, format, args);
  va_end (args);
  Info(0, str);
}

//with verbosity
void Message::Info(int level, const char *format, ...)
{
  if(level > _verbosity) return;
  char str[1024];
  va_list args;
  va_start (args, format);
  vsnprintf (str, 1024, format, args);
  va_end (args);
  fprintf(stdout, "Info    : %s\n", str);
  return;
}

//Warning...
void Message::Warning(const char *format, ...)
{
  char str[1024];
  va_list args;
  va_start (args, format);
  vsnprintf (str, 1024, format, args);
  va_end (args);
  Warning(0, str);
}

void Message::Warning(int level, const char *format, ...)
{
  if(level > _verbosity) return;
  char str[1024];
  va_list args;
  va_start (args, format);
  vsnprintf (str, 1024, format, args);
  va_end (args);
  //to write in bold (\33[1m) + red (\33[31m)
  const char *c0 = "", *c1 = "";
  c0 = "\33[1m\33[31m"; c1 = "\33[0m";
  //
  fprintf(stdout, "%sWarning : %s%s\n", c0,str,c1);
}

void Message::Check()
{
  Message::Info("Check param file...");
  Message::Info("-- Global params --");
  Message::Info("Res directory: %s", _resDir.c_str());
  Message::Info("Compute MC: %s", _ComputeMC?"Yes":"No");
  Message::Info("Post-Processing: %s", _Pos?"Yes":"No");
  Message::Info("-- Functions --");
  Message::Info("Number of function available: %d", _NFUN);
  for (int i =0; i<_NFUN; i++)
    {
      if(_FunChoice[i])  Message::Info("Function %d: Yes", i);
      else  Message::Info("Function %d: No (setting MC[%d]=0)", i);
      Message::Info("MC[%d] desired=%d %s", i, _desired_MC[i], _ComputeMC?"":"(Useless)");
    }
  Message::Info("-- Grid --");
  Message::Info("xi_min: %g", _xi_min);
  Message::Info("xi_max: %g", _xi_max);
  Message::Info("dxi: %g", _dxi);
  Message::Info("y_min: %g", _y_min);
  Message::Info("y_max: %g", _y_max);
  Message::Info("dy: %g", _dy);
}

//To quit properly like a boss
void Message::Finalize(int status)
{
#if defined(WITH_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif
  if(status == EXIT_SUCCESS)
    Message::Info("Exit with success");
  else
    Message::Warning("Exit with error (status %d)", status);
  exit(status);
}

//Parse param file to find the parameters wanted by the user ^_^
void Message::Parse()
{
  Message::Info("Parse param file \"%s\"...", _paramFile.c_str());
  std::ifstream pfile(_paramFile.c_str());
  if(!pfile.is_open())
    {
      Message::Warning("Paramameters file \"%s\" not found, exiting badly", _paramFile.c_str());
      Message::Finalize(EXIT_FAILURE);
    }
  else
    {
      std::string line;
      while (getline(pfile, line))
	{
	  //remove spaces
	  line.erase(remove_if(line.begin(), line.end(), isspace), line.end());
	  //remove comment
	  std::size_t found_comment = line.find("//");
	  if(found_comment >= 0 && found_comment < line.size()) line = line.substr(0, found_comment);
	  if(line.size() == 0) continue;
	  //check for equal sign
	  std::size_t found = line.find('=');
	  if(found > line.size()) Message::Warning("This line is unreadable: %s", line.c_str());
	  else
	    {
	      std::string keyword=line.substr(0,found);
	      std::string c_value=line.substr(found+1, line.size()-found);
	      //Check keyword
	      if(keyword == "resDir"){//set the "/" at the end of the folder directory's name (or not)
		if(c_value[c_value.size()-1] == '/') _resDir = c_value;
		else _resDir = c_value + "/";
	      }
	      //Check for functions, number of MC simulations...
	      for (int i=0; i<_NFUN; i++)
		{
		  std::ostringstream oss;
		  oss << i;
		  std::string mmc = "MC_" + oss.str();
		  std::string func = "FUN_" + oss.str();
		  int int_value = atoi(c_value.c_str());
		  if(keyword == mmc) { _desired_MC[i] = int_value;}
		  if(keyword == func){ _FunChoice[i] = (int_value ==0 ?0:1);}
		}
	      //Check for the grid !
	      if(keyword == "xi_min"){_xi_min = atof(c_value.c_str());}
	      if(keyword == "xi_max"){_xi_max = atof(c_value.c_str());}
	      if(keyword == "dxi"){ _dxi = atof(c_value.c_str());}
	      if(keyword == "y_min"){ _y_min = atof(c_value.c_str());}
	      if(keyword == "y_max"){ _y_max = atof(c_value.c_str());}
	      if(keyword == "dy"){ _dy = atof(c_value.c_str());}
	    }
	}
      pfile.close();
      //last check to put MC = 0 if function is not choiced (security ?)
      for(int i =0; i< _NFUN; i++)
	_desired_MC[i] = _desired_MC[i]*_FunChoice[i];
      if(_xi_min > _xi_max) Message::Warning("_xi_min > _xi_max, are you nuts ?");
      if(_y_min > _y_max) Message::Warning("_y_min > _y_max, are you nuts ?");
    }
}

void Message::DistributeWork(int nPointToDo, std::vector<int> *IndexOfPointToDo)
{
  //Given a number nPointToDo, it creates an array different for every MPI_Process. For example if there are 3 process...
  // Rank 0 : 0, 3, 6, 9, ...
  // Rank 1 : 1, 4, 7, 10, ...
  // Rank 2 : 2, 5, 8, 11, ...
  for (int i = _myRank; i < nPointToDo; i+=_nb_proc)
    IndexOfPointToDo->push_back(i);
}
