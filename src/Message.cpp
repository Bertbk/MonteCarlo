#include <iostream>
#include <fstream>
#include <sstream> //for osstream
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>
#include <time.h>
#include <math.h>
#include <algorithm> //for remove_if, isspace

#include <MonteCarlo/Message.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#if defined(HAVE_OMP)
#include <omp.h>
#endif


//using namespace std;
int Message::m_ComputeMC = 0;
int Message::m_Pos = 0;
int Message::m_Gmsh = 0;
int Message::m_verbosity = 4;
std::string Message::GMSHFileName = "res";
int Message::m_myRank = 0;
int Message::m_nb_proc = 1;
int Message::m_nb_threads = 1;
//PARAMETERS
//==============
const double Message::m_deuxpi = 8*atan(1.0);
double Message::m_lambda = 0.1;
//parameters of the oscillator
double Message::m_Y   = 1.0;
double Message::m_c0  = 1.0;
double Message::m_k   = 1.0;
//filtre
double Message::m_alpha = 0.5;  
//correlation of the noises between y and gamma (g)
double Message::m_ro = 0.1;
double Message::m_roc = sqrt(1.0 - m_ro*m_ro);
//parameters of the numerical procedure
double Message::m_T   = 500.0;
double Message::m_dt  = 0.0001;
double Message::m_sdt = sqrt(m_dt);
//Grid
double Message::m_xi_min=-5, Message::m_xi_max=5, Message::m_dxi = 5;
double Message::m_y_min=0, Message::m_y_max=5, Message::m_dy=5;
std::vector<double> Message::m_xi;
std::vector<double> Message::m_y;
//choice of function (for final computation)
const int Message::m_NFUN = 1;
std::vector<int> Message::m_FunChoice(m_NFUN);
std::vector<int> Message::m_desired_MC(m_NFUN);
std::string Message::m_paramFile = "param";
std::string Message::m_resDir = "res/";
std::string Message::m_helpDir = "help/";
int Message::m_restart = 0;
//File/folder names
std::string Message::DBext = ".db";
std::string Message::POSext = ".pos";
std::string Message::PointDatabase = "Points";
std::string Message::FullResRootName = "ResFun";
std::string Message::CurrentPointDatabase = "currentPoint";
std::string Message::PointFolderRootName = "Point";
std::string Message::FunResFolderRootName = "fun";
std::string Message::FunResRootName = "fun";
std::string Message::PointResRootName = "point_res_";
std::string Message::BackSlash = "/";
//Message
//-------
void Message::Initialize(int argc, char *argv[])
{
#if defined(HAVE_MPI)
  MPI_Init( &argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &m_myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &m_nb_proc);
  if(m_myRank == 0)
    Message::Info("Launched with MPI (%d processes)", m_nb_proc);
#endif
#if defined(HAVE_OMP)
#pragma omp parallel
  {
  m_nb_threads = omp_get_num_threads();
  }
  Message::Info("Launched with OpenMP (%d threads)", m_nb_threads);
#endif
  int i = 1;
  bool doCheckOnly = 0;
  bool showHelp = 1;
  while (i < argc) {
    if (argv[i][0] == '-') {
		if (!strcmp(argv[i] + 1, "check"))    { doCheckOnly = 1; i++; showHelp = 0; }
		else if (!strcmp(argv[i] + 1, "v"))   { m_verbosity = atof(argv[i + 1]); i += 2; showHelp = 0; }
		else if (!strcmp(argv[i] + 1, "par")) { m_paramFile = argv[i + 1]; i += 2; showHelp = 0; }
		else if (!strcmp(argv[i] + 1, "MC"))  { m_ComputeMC = 1; i++; showHelp = 0; }
		else if (!strcmp(argv[i] + 1, "pos")) { m_Pos = 1; i++; showHelp = 0; }
		else if (!strcmp(argv[i] + 1, "gmsh")) { m_Gmsh = 1; GMSHFileName = argv[i + 1]; i+=2; showHelp = 0; }
		else{ Warning("What the hell is this option (skipping) ? (%s)", argv[i] + 1); i++; }
	}
	else{ Warning("What the hell is this option (skipping) ? (%s)", argv[i]); i++; }
  }
  for(int i =0; i < m_NFUN; i++)
    {
      m_desired_MC[i] = 0;
      m_FunChoice[i] = 0;
    }
  if (showHelp){ Message::Help(); Message::Finalize(EXIT_SUCCESS); };
  //Parse param file
  Message::Parse();
  //Print info
  if(doCheckOnly) { Message::Check(); Message::Finalize(EXIT_SUCCESS);}
  else {if (m_verbosity > 0) Message::Check();}
  Message::BuildGrid();

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
  if(level > m_verbosity) return;
  char str[1024];
  va_list args;
  va_start (args, format);
  vsnprintf (str, 1024, format, args);
  va_end (args);
  fprintf(stdout, "Info[%d] : %s\n", m_myRank,str);
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
  if(level > m_verbosity) return;
  char str[1024];
  va_list args;
  va_start (args, format);
  vsnprintf (str, 1024, format, args);
  va_end (args);
  //to write in bold (\33[1m) + red (\33[31m)
  const char *c0 = "", *c1 = "";
  c0 = "\33[1m\33[31m"; c1 = "\33[0m";
  //
  fprintf(stdout, "%sWarning[%d]: %s%s\n", c0,m_myRank, str,c1);
}

//Debug...
void Message::Debug(const char *format, ...)
{
  char str[1024];
  va_list args;
  va_start (args, format);
  vsnprintf (str, 1024, format, args);
  va_end (args);
  Debug(0, str);
}

void Message::Debug(int level, const char *format, ...)
{
  if(level > m_verbosity) return;
  char str[1024];
  va_list args;
  va_start (args, format);
  vsnprintf (str, 1024, format, args);
  va_end (args);
  //to write in bold (\33[1m) + red (\33[31m)
  const char *c0 = "", *c1 = "";
  c0 = "\33[1m\33[34m"; c1 = "\33[0m";
  //
  fprintf(stdout, "%sDebug[%d]: %s%s\n", c0,m_myRank, str,c1);
}


// Show help of MonteCarlo (options, ...)
void Message::Help()
{
  if(m_myRank>0)
    return;
  std::cout << "Monte Carlo simulator\n";
  std::cout << "B. Thierry\n";
  std::cout << "Options: \n";
  std::cout << "  -par string         Select param file to parse (default = \"param\")\n";
  std::cout << "  -check              Check param file (nothing else is done)\n";
  std::cout << "  -v num              Set verbosity level (default = 4)\n";
  std::cout << "  -MC                 Launch the Monte Carlo computations\n";
  std::cout << "  -pos                Launch the Post Processing (can be used together with -MC or standalone)\n";
  std::cout << "  -gmsh  FILENAME     Print result on FILENAM.pos file (in the right folder) \n";
}

// Check the result of parsing the param file
void Message::Check()
{
  if(m_myRank !=0)
    return;
  Message::Info("Check param file...");
  Message::Info("========== Global params ==========");
  Message::Info("Res directory  : %s", m_resDir.c_str());
  Message::Info("Compute MC     : %s", m_ComputeMC?"Yes":"No");
  Message::Info("Post-Processing: %s", m_Pos?"Yes":"No");
  Message::Info("============ Functions ============");
  Message::Info("Number of function available: %d", m_NFUN);
  for (int i =0; i < m_NFUN; i++)
    {
      if(m_FunChoice[i])  Message::Info("Function %d: Yes", i);
      else  Message::Info("Function %d: No (setting MC[%d]=0)", i, i);
      Message::Info("MC[%d] desired=%d %s", i, m_desired_MC[i], m_ComputeMC?"":"(But not computation is asked)");
    }
  Message::Info("Restart: %d", m_restart);
  Message::Info("=============== Grid ==============");
  Message::Info("xi_min: %g", m_xi_min);
  Message::Info("xi_max: %g", m_xi_max);
  Message::Info("dxi   : %g", m_dxi);
  Message::Info("y_min : %g", m_y_min);
  Message::Info("y_max : %g", m_y_max);
  Message::Info("dy    : %g", m_dy);
  Message::Info("============ End Check ============");
}

//To quit properly like a boss
void Message::Finalize(int status)
{
#if defined(HAVE_MPI)
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
  if(Message::RootMpi())
    Message::Info("Parse param file \"%s\"...", m_paramFile.c_str());
  std::ifstream pfile(m_paramFile.c_str());
  if(!pfile.is_open())
    {
      Message::Warning("Paramameters file \"%s\" not found, exiting badly", m_paramFile.c_str());
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
		if(c_value[c_value.size()-1] == '/') m_resDir = c_value;
		else m_resDir = c_value + "/";
	      }
	      //Check for functions, number of MC simulations...
	      for (int i=0; i<m_NFUN; i++)
		{
		  std::ostringstream oss;
		  oss << i;
		  std::string mmc = "MC_" + oss.str();
		  std::string func = "FUN_" + oss.str();
		  int int_value = atoi(c_value.c_str());
		  if(keyword == mmc) { m_desired_MC[i] = int_value;}
		  if(keyword == func){ m_FunChoice[i] = (int_value ==0 ?0:1);}
		}
	      if(keyword == "RESTART"){m_restart = atoi(c_value.c_str());}
	      //Check for the grid !
	      if(keyword == "xi_min"){m_xi_min = atof(c_value.c_str());}
	      if(keyword == "xi_max"){m_xi_max = atof(c_value.c_str());}
	      if(keyword == "dxi"){ m_dxi = atof(c_value.c_str());}
	      if(keyword == "y_min"){ m_y_min = atof(c_value.c_str());}
	      if(keyword == "y_max"){ m_y_max = atof(c_value.c_str());}
	      if(keyword == "dy"){ m_dy = atof(c_value.c_str());}
	    }
	}
      pfile.close();
      //last check to put MC = 0 if function is not choiced (security ?)
      //???????????????????????????????????????????????????
      for(int i =0; i< m_NFUN; i++)
	m_desired_MC[i] = m_desired_MC[i]*m_FunChoice[i];
      if(m_xi_min > m_xi_max) Message::Warning("m_xi_min > m_xi_max, are you nuts ?");
      if(m_y_min > m_y_max) Message::Warning("m_y_min > m_y_max, are you nuts ?");
    }
}

void Message::BuildGrid()
{
  int nxi = ceil((m_xi_max - m_xi_min)/m_dxi) +1 ;
  int ny = ceil((m_y_max - m_y_min)/m_dy) + 1;
  int np = (nxi*ny);
  m_xi.reserve(np);
  m_y.reserve(np);
  for (int j = 0; j < ny ; j++)
    {
      for (int i = 0; i < nxi ; i++)
	{
	  m_xi.push_back(m_xi_min + i*m_dxi);
	  m_y.push_back(m_y_min + j*m_dy);
	}
    }
}

void Message::DistributeWork(int N, std::vector<int> *iStart, std::vector<int> *iEnd)
{
  //Given a number N, it creates two arrays :
  // iStart[myRank] = Starting index
  // iEnd[myRank] = Ending index (not comprise)
  // for (int i = iStart[myRank]; i < iEnd[myRank]; i++)
  //   ...
  iStart->resize(m_nb_proc);
  iEnd->resize(m_nb_proc);
  int step = N/m_nb_proc;
  int reste = N % m_nb_proc;
  (*iStart)[0] = 0;
  (*iEnd)[0] = step;
  if(reste > 0)
    (*iEnd)[0]++;
  for (int i = 1; i < m_nb_proc ; i++)
    {
      (*iStart)[i] = (*iEnd)[i-1];
      (*iEnd)[i] = (*iStart)[i] + step;
      if(i<reste)
	(*iEnd)[i] ++;
    }
}


int Message::GetThreadNum()
{
#ifdef HAVE_OMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

void Message::Barrier()
{
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  return;
}
