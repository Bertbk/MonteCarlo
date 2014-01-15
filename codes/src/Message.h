#ifndef _MESSAGE_H_
#define _MESSAGE_H_

#include <string>
#include "config.h"

// a class to manage messages
class Message{
 private:
  static int _ComputeMC; // Launch MC simulations ?
  static int _Pos; // Only (re-)compute average/std deviation ?
  static int _Gmsh;  // Print on Gmsh file
  static std::string _paramFile;
  static std::string _resDir;
  static int _verbosity;
  static int _myRank, _nb_proc;
  //---------------------------
  //Parameters of simulations
  //---------------------------
  //function (number available, which one are choosen, ...)
  static const int _NFUN; //this is hard coded !
  static std::vector<int> _FunChoice; //choice of the functions to compute final result
  //User wanted values
  static std::vector<int> _desired_MC;
  //constants
  static double _deuxpi;
  static double _lambda;
  //parameters of the oscillator
  static double _Y;
  static double _c0;
  static double _k;
  //filtre
  static double _alpha;
  //correlation of the noises between y and gamma (g)
  static double _ro;
  static double _roc;
  //parameters of the numerical procedure
  static double _T;
  static double _dt;
  static double _sdt;
  //Grid dimension
  static double _xi_min, _xi_max, _dxi, _y_min, _y_max, _dy;
 public:
  static void Initialize(int argc, char *argv[]);
  static void Info(int level, const char *format, ...);
  static void Info(const char *format, ...);
  static void Warning(const char *format, ...);
  static void Warning(int level, const char *format, ...);
  static int Precision(){return 17;} //set decimal precision for output file
  static int GetRank(){return _myRank;};
  static int GetNProc(){return _nb_proc;};
  static int GetComputeMC(){return _ComputeMC;};
  static int GetPos(){return _Pos;};
  static int GetGmsh(){return _Gmsh;};
  static int GetNFUN(){return _NFUN;};
  static std::string GetResDir(){return _resDir;};
  static void Check();
  static void Parse();
  static void Finalize(int status);
  //=======================
  // PARAMETERS FUNCTIONS
  //=======================
  static double GetDeuxPi(){return _deuxpi;};
  static double GetLambda(){return _lambda;};
  static double GetY(){return _Y;};
  static double GetC(){return _c0;};
  static double GetK(){return _k;};
  static double GetAlpha(){return _alpha;};
  static double GetRo(){return _ro;};
  static double GetRoc(){return _roc;};
  static double GetT(){return _T;};
  static double GetDt(){return _dt;};
  static double GetSdt(){return _sdt;};
  //MPI DISTRIBUTER
  static void DistributeWork(int nPointToDo, std::vector<int> *IndexOfPointToDo);
};

#endif
