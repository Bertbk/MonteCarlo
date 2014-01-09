#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <time.h>
#include <stdarg.h>
#include <math.h>
#include <string>
#include <string.h>
//#include <cstring>

#include "Message.h"

#if defined(WITH_MPI)
#include "mpi.h"
#endif

#if defined(WITH_OMP)
#include "omp.h"
#endif


using namespace std;

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

//choice of function (for final computation)
int  Message::NFUN = 1;
std::vector<int> Message::_FunChoice(NFUN);
std::string Message::_paramFile = "param";
std::string Message::_resDir = "res";

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
  while (i < argc) {
    if (argv[i][0] == '-') {
      if (!strcmp(argv[i]+1, "v"))          { _verbosity = atof(argv[i+1]) ; i+=2 ; }
      else if (!strcmp(argv[i]+1, "check")) { Message::Check(); Message::Finalize(EXIT_SUCCESS);}
      else if (!strcmp(argv[i]+1, "par"))   { _paramFile = argv[i+1]; i+=2; }
    }
    else
      {
	printf("What the hell is this option ?\n");
	i++;
      }
  }
  //Parse param file
  Message::Parse();
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
  Message::Info("Check...");  
}

void Message::Parse()
{
  Message::Info("Parse param file \"%s\"...", _paramFile.c_str());
}


void Message::Finalize(int status)
{
#if defined(WITH_MPI)
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif
  if(status == EXIT_SUCCESS)
    Message::Info("Exit with success");
  else
    Message::Info("Exit with error (status %d)", status);
  exit(status);
}
