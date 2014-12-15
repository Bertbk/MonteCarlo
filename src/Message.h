#ifndef _MESSAGE_H_
#define _MESSAGE_H_

#include <string>


// a class to manage messages
class Message{
 private:
  static int m_ComputeMC; // Launch MC simulations ?
  static int m_Pos; // Only (re-)compute average/std deviation ?
  static int m_Gmsh;  // Print on Gmsh file
  static std::string m_paramFile;
  static std::string m_resDir;
  static int m_verbosity;
  static int m_myRank, m_nb_proc;
  //---------------------------
  //Parameters of simulations
  //---------------------------
  //function (number available, which one are choosen, ...)
  static const int m_NFUN; //this is hard coded !
  static std::vector<int> m_FunChoice; //choice of the functions to compute final result
  //User wanted values
  static std::vector<int> m_desired_MC;
  //constants
  static double m_deuxpi;
  static double m_lambda;
  //parameters of the oscillator
  static double m_Y;
  static double m_c0;
  static double m_k;
  //filtre
  static double m_alpha;
  //correlation of the noises between y and gamma (g)
  static double m_ro;
  static double m_roc;
  //parameters of the numerical procedure
  static double m_T;
  static double m_dt;
  static double m_sdt;
  //Grid dimension
  static double m_xi_min, m_xi_max, m_dxi, m_y_min, m_y_max, m_dy;
 public:
  static void Initialize(int argc, char *argv[]);
  static void Info(int level, const char *format, ...);
  static void Info(const char *format, ...);
  static void Warning(const char *format, ...);
  static void Warning(int level, const char *format, ...);
  static int Precision(){return 17;} //set decimal precision for output file
  static int GetRank(){return m_myRank;};
  static int GetNProc(){return m_nb_proc;};
  static int GetComputeMC(){return m_ComputeMC;};
  static int GetPos(){return m_Pos;};
  static int GetGmsh(){return m_Gmsh;};
  static const int GetNFUN(){return m_NFUN;};
  static std::string GetResDir(){return m_resDir;};
  static void Check();
  static void Parse();
  static void Finalize(int status);
  static void Help();
  //=======================
  // PARAMETERS FUNCTIONS
  //=======================
  static double GetDeuxPi(){return m_deuxpi;};
  static double GetLambda(){return m_lambda;};
  static double GetY(){return m_Y;};
  static double GetC(){return m_c0;};
  static double GetK(){return m_k;};
  static double GetAlpha(){return m_alpha;};
  static double GetRo(){return m_ro;};
  static double GetRoc(){return m_roc;};
  static double GetT(){return m_T;};
  static double GetDt(){return m_dt;};
  static double GetSdt(){return m_sdt;};
  //MPI DISTRIBUTER
  static void DistributeWork(int nPointToDo, std::vector<int> *IndexOfPointToDo);
};

#endif
