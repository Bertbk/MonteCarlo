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
  std::vector <Point*> Points;
  std::vector <Point*> PointsToDo;
  static std::string DBext;
  static std::string PointDatabase;
  static std::string FullResRootName;
  static std::string CurrentPointDatabase;
  static std::string FunResFolderRootName;
  static std::string FunResRootName;
  static std::string PointFolderRootName;
  static std::string PointResRootName;
  int Npoints, NpointsToDo, MaxId;

  //Parse the Root files: points.db, resfunXX.db, ...
  void ParseRootFiles();
  //loop on PointXX folder and parse the files and funXX folder...
  void ParsePointFiles();

  void BuildFolderPoint(int id); // Create a folder for the point with empty files
  void RebuildPointsDb(); // Rebuild points.db
 public:
  //Constructor
  Database(std::string resdir);
  ~Database();
  //Init Database: read files, build PointsDone vector, ...
  void Init();
  //Provide an array of pointers to the Points to be treated by the MC solver
  void UpdatePointsToDo(std::vector<double> X, std::vector<double> Y, std::vector<int> MC_To_Do);
  int FindPoint(double x, double y);
  Point* GetPointToDo(int index);
  //Print information on every points
  void PrintPoints();

};

#endif
