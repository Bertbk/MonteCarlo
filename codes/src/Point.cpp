#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<stdio.h>
#include <vector>
#include<time.h>

#include "Message.h"
#include "Point.h"

using namespace std;

void Point::ReadAllPoints(std::vector<Point*> *PointDone)
{
  //Read all folder and build the right Point
}

void Point::CreatePointsToDo(std::vector<Point*> *PointToDo, std::vector<Point*> *PointDone)
{
  //Use the param file to build point according to :
  // -- The given grid (param)
  // -- The already existing point (PointDone)
  //Destruction of the point must be achieved by the user !
}


/*

//constructor
Point::Point()
{
  _id = -1;
  _xi = -1.;
  _y = -1.;
  _MC = -1;
  _npoints ++;
}

Point::Point(double xi, double y, int MC){
  _id = -1;    
  _xi = xi;
  _y = y;
  _MC = MC;
  _npoints ++;
}

Point::Point(int id, double xi, double y, int MC){
  if(id > 0)
    _id = id;
  else
    _id = -1;    
  _xi = xi;
  _y = y;
  _MC = MC;
  _npoints ++;
}

//Change parameters
void Point::Set(int id ,double xi, double y, int MC)
{
  if(_id != id)
    {
      if(_id != -1)
	Message::Warning("[Point id %d] Change of id (?!) from %g to %g !", _id, _id, id);
      _id = id;
    }
  Set(xi, y, MC);
}
//Change parameters
void Point::Set(double xi, double y, int MC)
{
  if(_xi != xi)
    {
      Message::Warning("[Id %d] xi-coordinate changes from %g to %g !", _id, _xi, xi);
      _xi = xi;
    }
  if(_y != y)
    {
      Message::Warning("[Id %d] y-coordinate changes from %g to %g !", _id, _y, y);
      _y = y;
    }
  if(_MC != MC)
    {
      Message::Warning("[Id %d] MC changes from %d to %d !", _id, _MC, MC);
      _MC = MC;
    }
  return;
}

//Check if this point has already been done or not
void Point::Check()
{
  int nres = MyResults::GetNRES();
  for (int i = 0; i < nres; i++)
    {
      double x=MyResults::GetX(i), y = MyResults::GetY(i);
      if(_xi == x && _y == y)
      {
	int id = MyResults::GetId(i), MC = MyResults::GetMC(i);
	int MC_to_do = _MC - MC;
	this->SetId(id);
	if(MC_to_do <= 0)
	  Message::Warning("[Id %d] Already exist with desired number of simulations", id);
	else
	  Message::Warning("[Id %d] Already exist, only %d simu to reach %d", id, MC_to_do, _MC);
	this->Set(_xi, _y, MC_to_do);
      }
    }  
  if(_id == -1) //Set new id
    this->SetNewId();
  return;
}

//Ask MyResults for a new id and set it to the current point
void Point::SetNewId()
{
  int id = MyResults::GetNewId();
  _id = id;
  Message::Info(2,"New ID provided");
  return;
}


void Point::PrepareMyFile(char *res_dir, int *file_id)
{
  std::cout.precision(Message::Precision());
  int id = _id;
  double xi=_xi, y=_y;
  char command[128], res_dir_loc[128];
  sprintf(res_dir_loc, "./%s/Id%d/", res_dir, id);
  sprintf(command, "if ! test -d %s; then mkdir %s; fi", res_dir_loc, res_dir_loc);
  //  sprintf(command, "mkdir %s", res_dir_loc);
  system(command);
  //create file (search file id that is free (no erase))
  int fcpt = 0;
  char file_aux[128];
  sprintf(file_aux, "%sres%d.mc", res_dir_loc, fcpt);
  std::ifstream fNewRes_aux(file_aux);
  while(!fNewRes_aux.fail())
    {
      fNewRes_aux.close();
      fcpt++;
      sprintf(file_aux, "%sres%d.mc", res_dir_loc, fcpt);
      fNewRes_aux.open(file_aux);
    }
  *file_id = fcpt;
  
  //create/verify file containing (x,y) coordinate
  char coord_file[128];
  sprintf(coord_file, "%smyCoord.mc", res_dir_loc);
  //small checking on the coordinate of the points...
  std::ifstream fcoord_in(coord_file);
  if((fcoord_in.fail()))
    {
      std::ofstream fcoord(coord_file);
      fcoord << xi << "\n" << y;
      fcoord.close();
    }
  else{//check file
    double xi_aux, y_aux;
    fcoord_in >> xi_aux >> y_aux;
    fcoord_in.close();
    bool test_file=(xi == xi_aux && y == y_aux);
    if(!test_file)
      {//There is a (huge) problem, better keep all information in a separate file
	Message::Warning("[Id %d] Bad coordinate detected !", id);
	char coord_file_aux[128];
	sprintf(coord_file_aux, "%s_aux%d", coord_file, fcpt);
	std::ofstream fcoord(coord_file_aux);
	fcoord << xi << "\n" << y;
	fcoord.close();
      }
  }
  return;
}

void Point::LaunchMC(char *traj_dir, int seed)
{
  Message::Info("[Proc %d] I will do %d MC tests on point %g %g", Message::GetCommRank(), _MC, _xi, _y);

  /*
  //change seed
  //time(NULL) - 360000*myRank
  srand(seed);
  //Create directory (if doesn't exist)
  int file_id;
  PrepareMyFile(traj_dir, &file_id);
  char traj_dir_loc[128], traj_file[128];
  sprintf(traj_dir_loc, "./%s/Id%d/", traj_dir, _id);
  sprintf(traj_file, "%sres%d.mc", traj_dir_loc, file_id);
  std::ofstream fNewRes(traj_file);
  
  //#pragma omp parallel for private(imc)
  myStr << "$MC\n"<< _MC << "\n";
  for (int imc = 0 ; imc < _MC ; imc++)
    {
      Message::Info("%d",imc);
      std::vector<std::vector<double> > traj_imc;
      ShortCyclePlus(&traj_imc);
      std::stringstream myStr(ostringstream::out|ios::binary);
      myStr << traj_imc.size() << "\n";
      for (int i = 0; i < traj_imc.size(); i++)
	myStr << std::setprecision(Message::Precision()) << traj_imc[i][0] << " " << std::setprecision(Message::Precision()) << traj_imc[i][1] << "\n";      
      fNewRes << myStr.str();
    }
  fNewRes.close(); 
  */
/*}



void Point::ShortCyclePlus(std::vector<std::vector<double> > *trajectories)
{
  //run a trajectory starting from (xi,y) \in DELTAPLUS of the solution on the boundary.
  double integ = 0.0;
  int cpt = 0;
  double t = 0.;
  double T = Message::GetT();
  double dt = Message::GetDt();
  double sdt = Message::GetSdt();
  double ro = Message::GetRo();
  double roc = Message::GetRoc();
  double c0 = Message::GetC();
  double k = Message::GetK();
  double Y = Message::GetY();
  double alpha = Message::GetAlpha();
  int it_max = T/dt;
  double xi = _xi, y = _y;
  (*trajectories).resize(it_max);
  for(int it = 0 ; it < it_max ; it++)
    {
      (*trajectories)[it].resize(2);
      t += dt;
      if(t > T){ t = T;}
      //noise
      double g1 = gauss();
      double g2 = gauss();
      double xi_aux = xi, y_aux = y;
      xi += -alpha*xi_aux*dt + sdt*g1;
      y  += -(alpha*xi_aux + c0*y_aux + k*Y)*dt + sdt*(ro*g1 + roc*g2);
      //Save value of (xi,y) for further computation
      (*trajectories)[cpt][0]=xi;
      (*trajectories)[cpt][1]=y;
      cpt++;
      if (y<0) 
	break;
    }
  (*trajectories).resize(cpt);
  return;
}
*/
