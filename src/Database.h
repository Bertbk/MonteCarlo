#ifndef _Database_H_
#define _Database_H_

#include <vector>
#include <string>
#include <math.h>
#include <algorithm>
#include "Message.h"
#include "Point.h"

class Database {
 private:
  std::string m_resDir;
  std::vector <std::vector < int > > NResByFun; // Number of results per function (per points)
  std::vector <std::vector < double > > MeanByFun;// Average of the results (mean) per function (per points)
  std::vector <std::vector < double > > StdDevByFun;// Standard deviation of the results per function (per points)
  std::vector <Point*> PointsDone;
  std::vector <Point*> PointsToDo;
  static std::string DBext;
  static std::string PointDatabase;
  static std::string FullResRootName;
  static std::string CurrentPointDatabase;
  static std::string FunResFolderRootName;
  static std::string FunResRootName;
  static std::string PointFolderRootName;
  static std::string PointResRootName;
  int NpointsDone, NpointsToDo, MaxId;

  //Parse the Root files: points.db, resfunXX.db, ...
  void ParseRootFiles();
  //loop on PointXX folder and parse the files and funXX folder...
  void ParsePointFiles();

 public:
  //Constructor
  Database(std::string resdir);
  ~Database();
  //Init Database: read files, build PointsDone vector, ...
  void Init();

  Point* GetPointToDo(int index);
  //Print information on every points
  void PrintPointsDone();

};

#endif
