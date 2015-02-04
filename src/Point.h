#ifndef _Point_H_
#define _Point_H_

#include<math.h>
//#include<vector>

//Class that contains every points to be computed 
class Point{
 private:
  int m_id; //id of the point
  std::string m_id_str; // id, version string
  std::string m_IdDir;
  std::string m_myDir;
  double m_xi, m_y;
  std::vector<int> m_MC; //number of simu already done (for each function)
  std::vector<int> m_MC_to_do; //number of simu to do (for each function)
  std::vector<int> m_NResFiles; //number of files already written (for each function)
  std::vector<double> m_average, m_stddev; //current result in memory

  static std::string BackSlash;
  static std::string DBext;
  static std::string PointDatabase;
  static std::string FullResRootName;
  static std::string CurrentPointDatabase;
  static std::string FunResFolderRootName;
  static std::string FunResRootName;
  static std::string PointFolderRootName;
  static std::string PointResRootName;

  void WriteOnFile(std::vector<std::vector<double>*> *results);

 public:
  Point(int id, double xi, double y);
  Point(Point *p);
  static void ReadAllPoints(std::vector<Point*> *PointDone);
  static void CreatePointsToDo(std::vector<Point*> *PointToDo, std::vector<Point*> *PointDone);
  // Set some values
  void SetMCDone(int ifun, int MC){m_MC[ifun] = MC;}
  void SetMC_to_do(int ifun, int MC_to_do){m_MC_to_do[ifun] = MC_to_do;}
  void SetNResFiles(int ifun, int nfiles){m_NResFiles[ifun] = nfiles;}
  void SetAverage(int ifun, int average){m_average[ifun] = average;}
  void SetStdDev(int ifun, int stddev){m_stddev[ifun] = stddev;}
  //Monte Carlo simulations
  void SetMCToDo(std::vector<int> *desired_MC); // Set the right MC_to_do, by comparing (differencing) m_MC and m_Desired_MC...
  void LaunchMC();
  void ShortCyclePlus(std::vector<double> *integrals);
  //double uniform(){return (1+2)/(1+(double)RAND_MAX);   }
  static double uniform(){return (1+(double)rand())/(1+(double)RAND_MAX);   }
  static double gauss(){return sqrt(-2.*log(uniform()))*cos(Message::GetDeuxPi()*uniform());}
  static double f(double xi, int i);
  static double gplus(double xi, double y, int i);
  void Print();
  //Get functions
  int GetId(){return m_id;}
  double GetXi(){return m_xi;}
  double GetY(){return m_y;}
  int GetMC(int ifun){return m_MC[ifun];}
  int GetMCToDo(int ifun){return m_MC_to_do[ifun];}
  double GetAverage(int ifun){return m_average[ifun];}
  double GetStdDev(int ifun){return m_stddev[ifun];}
  int GetNResFiles(int ifun){return m_NResFiles[ifun];}
  std::string GetIdDir(){return m_IdDir;}
  /*  static int m_npoints;
  //=========================
  //TRAJECTORIES COMPUTATION
  //=========================
  //trajectories function
  static double uniform(){return (1+(double)rand())/(1+(double)RAND_MAX);};
  //  static double uniform(){return (1+2)/(1+(double)RAND_MAX);};
  static double gauss(){return sqrt(-2.*log(uniform()))*cos(Message::GetDeuxPi()*uniform());};
  void ShortCyclePlus(std::vector<std::vector<double> > *trajectories);
  //launch the monte carlo simulations
  void LaunchMC(char *traj_dir, int seed);
  /*
  //=====
  //DATA
  //=====
  static double f(double xi){ return exp(-xi*xi);};
  static double gplus(double xi, double y){return f(xi)*f(y);};
  */
};

#endif
