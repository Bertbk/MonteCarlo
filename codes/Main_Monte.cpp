// 2D Monte Carlo Method
// Results are stored in a new directory
// and are printed in two different matlab files
// for i=0,..., Nx
//   for j=0,..., Ny
//     res_i_?_j_?.m : the MC "final" results for the given functions f, g, ...
//     traj_i_?_j_?.m ; the MC trajectories. Usefull to compute final results by changing f and g without having to relaunch Monte Carlo !

#include<iostream>
#include<cstdlib>
#include<ctime>
#include<cmath>
#include<fstream>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>
#include <vector>
#include <iomanip>
#include <sstream>


/*#include<mpi.h>
  #include<omp.h>*/

#include "Message.h"
#include "MyResults.h"
#include "Point.h"

using namespace std;


// function declaration
//double shortcycleplus(double xi, double y, double lambda, std::vector<std::vector<double> > *res);

//probabilistic variables

/*const double deuxpi = 8*atan(1.0);
double uniform(){return (1+(double)rand())/(1+(double)RAND_MAX);   }
//double uniform(){return (1+2)/(1+(double)RAND_MAX);   }
double gauss(){return sqrt(-2.*log(uniform()))*cos(deuxpi*uniform());}
//DATA
double f(double xi){ return exp(-xi*xi);}
double gplus(double xi, double y){return f(xi)*f(y);}
*/

/*//parameters of the oscillator
const double Y   = 1.0;
const double c0  = 1.0;
const double k   = 1.0;

//filtre
const double alpha = 0.5;

//correlation of the noises between y and gamma (g)
const double ro = 0.1;
const double roc = sqrt(1.0-ro*ro);

//parameters of the numerical procedure
const double T   = 500.0;
const double dt  = 0.0001;
const double sdt = sqrt(dt);
*/
int main(int argc, char *argv[])
{
  //initiatilization (reading arguments, launching MPI,...)
  Message::Initialize(argc, argv);
  int myRank=Message::GetCommRank();
  int nb_proc=Message::GetNProc();
  std::cout.precision(Message::Precision());
  char res_dir[128], command[128];
  sprintf(res_dir, "res");
  //create res dir (if doesn't exist)
  //  sprintf(command, "mkdir %s", res_dir);
  sprintf(command, "if ! test -d %s; then mkdir %s; fi", res_dir, res_dir);
  system(command);

  //----------------------------------------------------------------
  //GRID POINTS TO COMPUTE
  //  double lambda = 0.1;
  //Number of Monte Carlo simulations that we want to do
  int n_MC = 5;

  //size of the domain
  double xi_min = -5.0; double xi_max = 5.0;
  double y_min = 0;     double y_max = 5.0;
  double Lxi = xi_max - xi_min;
  double Ly = y_max - y_min;

  //Number of point in each direction
  int Nxi = 2;
  int Ny = 2;
  int Ntot = Nxi*Ny; //total number of points

  //discretization step
  double hxi = Lxi/Nxi;
  double hy = Ly/Ny;

  //table of coordinates
  std::vector<double> tab_xi, tab_y;
  tab_xi.resize(Nxi); tab_y.resize(Ny);
  for(int ixi = 0; ixi < Nxi ; ixi ++) tab_xi[ixi] = (ixi + 0.5)*hxi;
  for(int iy = 0; iy < Ny ; iy ++) tab_y[iy] = (iy + 0.5)*hy;
  Message::Info("Nxi = %d, Ny = %d.", Nxi, Ny);
  //----------------------------------------------------------------
  //Reading old results
  MyResults::Initialize(argc, argv);
  //Create points and distribute them to MPI Process
  std::vector<Point> MyPoints(Nxi*Ny);
  int n_points = 0; // number of point I will do (local)
  int cpt = 0;
  for(int ixi = 0; ixi < Nxi; ixi++)
    {
      for(int iy = 0; iy < Ny; iy++)
	{
	  Point a(tab_xi[ixi], tab_y[iy], n_MC);
	  a.Check(); //has it already been done ? if yes, it modifies this->MC in consequence
	  if(a.GetMC()>0)
	    {
	      //Some computation must be done for this point to reach MC number of simus
	      if(cpt % nb_proc == myRank)
		{
		  MyPoints[n_points] = a; //...and it is for me !
		  n_points ++;
		}
	      cpt++;
	    }
	}
    }
  MyPoints.resize(n_points);
  Message::Info("[Proc %d] I will do %d points", myRank, n_points);
  // Now the point are distributed throughout the MPI Process.

  //What we need is to solve them ...
  //#pragma omp parallel
    // {
    ////little check
    //if(omp_get_thread_num() == 0)
  // printf("Proc %d have %d threads\n", myRank, omp_get_num_threads());
  //}

  //ensure that directory has been created and everybody is ready to work
  //  MPI_Barrier(MPI_COMM_WORLD);
  clock_t tstart = clock();
  clock_t tstarts = time(NULL);

  for (int iP = 0; iP < n_points; iP++)
    //    {
    MyPoints[iP].LaunchMC(res_dir, time(NULL) - 360000*myRank);


      /*      Point a = MyPoints[iP];
      //change seed
      srand(time(NULL) - 360000*myRank);
      double xi, y;
      int MC_loc, id;
      a.GetValues(&id, &xi, &y, &MC_loc);
      //      Message::Info("[Proc %d] I will do %d MC tests on point %g %g", myRank, MC_loc, xi, y);
      //Create directory (if doesn't exist)
      int file_id;
      a.PrepareMyFile(res_dir, &file_id);
      char res_dir_loc[128], res_file[128];
      sprintf(res_dir_loc, "./%s/Id%d/", res_dir, id);
      sprintf(res_file, "%sres%d.mc", res_dir_loc, file_id);
      std::ofstream fNewRes(res_file);

      //omp_lock_t writelock;
      //OMP_INIT_LOCK(&writelock);

      //#pragma omp parallel for private(imc) reduction(+:sum_vec) reduction(+:square_sum_vec)
      for (int imc = 0 ; imc < MC_loc ; imc++)
	{
	  Message::Info("%d",imc);
	  std::vector<std::vector<double> > res_imc;
	  shortcycleplus(xi, y, lambda, &res_imc);
	  std::stringstream myStr(ostringstream::out|ios::binary);
	  myStr << "$MC\n"<< MC_loc << "\n";
	  myStr << res_imc.size() << "\n";
	  for (int i = 0; i < res_imc.size(); i++)
	    myStr << std::setprecision(Message::Precision()) << res_imc[i][0] << " " << std::setprecision(Message::Precision()) << res_imc[i][1] << "\n";

	  //OMP_SET_LOCK(&writelock);
	  fNewRes << myStr.str();
	  //OMP_UNSET_LOCK(&writelock);
	}
      //OMP_DESTROY_LOCK(&writelock);
      fNewRes.close();
    }
      */
  //MERGE RESULTS FILE
  //End HAPPILY.
  //  MPI_Finalize();
  Message::Info("Program ended");
  return 0;

  /*  
  //Divide the work between MPI process
  int ind_start, ind_end; //start/end indices
  int Ntot_loc = Ntot/nb_proc; // local number of points (=jobs)
  int N_rest = Ntot%nb_proc; // non assigned points
  if(N_rest !=0 && myRank < N_rest)
    Ntot_loc ++; // this process will work one time more than some others
  ind_start = myRank*Ntot_loc;
  if(myRank >= N_rest)//previous processus has one more job to do, I must start "later"
    ind_start += min(myRank, N_rest);
  ind_end = ind_start + Ntot_loc;
  printf("Hi, Proc %d here : Ntot = %d, Ntot_loc = %d, ind_start = %d and ind_end = %d\n", myRank, Ntot, Ntot_loc, ind_start, ind_end);

#pragma omp parallel
{
  //little check
  if(omp_get_thread_num() == 0)
    printf("Proc %d have %d threads\n", myRank, omp_get_num_threads());
}
  //number of data to be stored in res_ files
  int n_data = 5; //X,Y,Z, Moments of order 1 and 2
  //vector of results
  std::vector<std::vector<double> > value_loc(Ntot_loc);
  int cpt_value = 0; //counter in value_loc
  //ensure that directory has been created and everybody is ready to work
  MPI_Barrier(MPI_COMM_WORLD);
  clock_t tstart = clock();
  clock_t tstarts = time(NULL);
  for (int i_loc = ind_start ; i_loc < ind_end ; i_loc++)
    {
      //One different seed per point and MPI process to "ensure" randomness
      srand(time(NULL) - 360000*myRank);
      //global index of current point
      int ixi, iy, iz;
      LocalToGlobal(i_loc, Ny, Nz, &ixi, &iy, &iz);
      //coordinates of the current point (crt)
      */      /*
      // WHY O.5 ??!!!
      double xi_i = -Lxi + (ixi + 0.5)*hg;
      double yj = (iy + 0.5)*hy;
      double zk = 0; //
	      *//*
      double xi_crt = tab_xi[ixi]; //xi_min + (ixi + 0.5)*hg;
      double y_crt = tab_y[iy]; //y_min + (iy + 0.5)*hy;
      double z_crt = tab_y[iz]; //z_min + iz*hz; //weirdo

      printf("Proc %d doing \ti_loc = %d, \tix = %d, \tiy = %d, \txi_i = %g, \tyj = %g\n", myRank, i_loc, ixi, iy, xi, y);
      //create files containing all the result for that point and the trajectories
      FILE *fid_loc;
      char filename[50];
      //RESULT file
      sprintf(filename, "res/res_i_%d_j_%d_k_%d.m", ix, iy, iz);
      fid_res = fopen(filename, "w");
      fprintf(fid_res, "\%\% RESULT FOR A TABLE OF %d x %d x %d", Nxi, Ny, Nz);
      fprintf(fid_res, "%%%% ixi = %d iy = %d iz = %d\n", ixi, iy, iz);
      fprintf(fid_res, "xi = 0.16%g\n", xi_crt);
      fprintf(fid_res, "y = 0.16%g\n", y_crt);
      fprintf(fid_res, "z = 0.16%g\n", z_crt);
      fprintf(fid_res, "res = [\n");
      //TRAJECTORY file
      sprintf(filename, "res/traj_i_%d_j_%d_k_%d.m", ix, iy, iz);
      fid_traj = fopen(filename, "w");
      fprintf(fid_res, "\%\% RESULT FOR A TABLE OF %d x %d x %d", Nxi, Ny, Nz);
      fprintf(fid_traj, "%%%% ixi = %d iy = %d iz = %d\n", ixi, iy, iz);
      fprintf(fid_traj, "xi = 0.16%g\n", xi_crt);
      fprintf(fid_traj, "y = 0.16%g\n", y_crt);
      fprintf(fid_traj, "z = 0.16%g\n", z_crt);
      fprintf(fid_traj, "traj = [\n");

      double sum_vec = 0, square_sum_vec = 0;
      int imc;
      //MC tests done on multiple threads with shared memory (OpenMP)
#pragma omp parallel for private(imc) reduction(+:sum_vec) reduction(+:square_sum_vec)
      for (imc = 0 ; imc < MC ; imc++)
	{
	  double res_imc = shortcycleplus(xi_crt, y_crt, lambda);
	  fprintf(fid_res, "%.16g\n", res_imc);
	  sum_vec = sum_vec + res_imc; //reduction
	  square_sum_vec = square_sum_vec + res_imc*res_imc; //reduction
	}
      //close file with MC results for this point
      fprintf(fid_res, "];\n");
      fclose(fid_loc);
      //Compute value of interest
      double X_1 = sum_vec/MC, X_2 = square_sum_vec/MC;
      value_loc[cpt_value].resize(n_data);
      value_loc[cpt_value][0] = xi_i;
      value_loc[cpt_value][1] = yj;
      value_loc[cpt_value][2] = zk;
      value_loc[cpt_value][3] = X_1; //average
      value_loc[cpt_value][4] = X_2; // moment order 2
      cpt_value++;
    }
  //Everything is done, now a file summarizing every result must be writen
  double timetot = difftime(clock(), tstart)/CLOCKS_PER_SEC;
  double timetots = difftime(time(NULL), tstarts);
  printf("Proc %d finishes its %d jobs in %g seconds (cpu %g)!\n", myRank,  Ntot_loc, timetots, timetot);
		*/  /*
  // Every process will send their results "value_loc" to master proc
  std::vector<int> tab_Ntot_loc(nb_proc);
  MPI_Gather(&Ntot_loc, 1, MPI_INT, &tab_Ntot_loc[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
  if(myRank > 0) //send
    {
      for (int ipoint = 0 ; ipoint < Ntot_loc ; ipoint ++)
	MPI_Send(&value_loc[ipoint][0], n_data, MPI_DOUBLE, 0, ipoint, MPI_COMM_WORLD);
    }else //rank == 0: receive everything 
    {
      //file containing all the results
      FILE *fid_res;
      fid_res = fopen("res/MATLAB_res.m", "w");
      fprintf(fid_res, "\%\% RESULT FOR A TABLE OF %d x %d x %d", Nxi, Ny, Nz);
      fprintf(fid_res, "\nres_monte = [\n");
      //receive
      for (int irank = 0 ; irank < nb_proc ; irank ++)
	{
	  for (int ipoint = 0 ; ipoint < tab_Ntot_loc[irank] ; ipoint ++)
	    {
	      std::vector<double> ValueAux(n_data);
	      MPI_Status status;
	      if(irank == 0)
		ValueAux = value_loc[ipoint];
	      else
		MPI_Recv(&ValueAux[0], n_data, MPI_DOUBLE, irank, ipoint, MPI_COMM_WORLD, &status);
	      //write values on file
	      for(int idata = 0 ; idata < n_data ; idata ++)
		fprintf(fid_res,"%.16g ", (double)ValueAux[idata]);
	      fprintf(fid_res,";\n");
	    }
	}
      fprintf(fid_res,"];\n");
      fclose(fid_res);
    }
		    *//*
  //ending
  if(myRank == 0)
    printf("Everything is finished\n");
  MPI_Finalize();*/
  return 0;
}

/*
//
double shortcycleplus(double xi, double y, double lambda, std::vector<std::vector<double> > *res)//z=Y,y>0
{
  //run a trajectory starting from (xi,y) \in DELTAPLUS of the solution on the boundary.
  //PLUS save the trajectory in a file
  /////////////////////////////////////////////////////////////////////////////////////
  double integ = 0.0;
  int cpt = 0;
  double t = 0.;
  int it_max = T/dt;
  (*res).resize(it_max);
  //  for(double t = dt ; t <=T ; t+=dt)
  for(int it = 0 ; it < it_max ; it++)
    {
      (*res)[it].resize(2);
      t +=dt;
      if(t>T){t=T;}
      //noise
      double g1 = gauss();
      double g2 = gauss();
      double xi_aux = xi, y_aux = y;
      xi += -alpha*xi_aux*dt + sdt*g1;
      y  += -(alpha*xi_aux + c0*y_aux + k*Y)*dt + sdt*(ro*g1 + roc*g2);
      //Save value of (xi,y) for further computation
      (*res)[cpt][0]=xi;
      (*res)[cpt][1]=y;
      cpt++;
      ///////////////////////////////////////
      integ += exp(-lambda*t)*gplus(xi,y);
      if (y<0) 
	break;
    }
  (*res).resize(cpt);
  return exp(-lambda*t)*f(xi) + integ*dt;
//  return -10;
}

*/
