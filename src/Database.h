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
  std::vector <Point*> PointsDone;
  std::vector <Point*> PointsToDo;
  static std::string DBext;
  static std::string PointDatabase;
  static std::string FullResRootName;
  static std::string CurrentPointDatabase;
  static std::string FunResRootName;
  static std::string PointResRootName;
  int NpointsDone, NpointsToDo;
 public:
  //Constructor
  Database(std::string resdir);
  ~Database();
  //Init Database: read files, build PointsDone vector, ...
  void Init();

  Point* GetPointToDo(int index);
  void ParseRootFiles();
  void PrintPointsDone();

};

#endif
