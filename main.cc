#include <iostream>
#include <sstream>
#include <mpi.h>
#include "parallel.h"

#include "partitioning.h"
#include "grid.h"
#include "output_vtk.h"
#include "linear_system.h"
#include "conjugate_gradient.h"
#include "assembling.h"

using namespace std;

int main( int argc, char** argv)
{
	MPI_Init(&argc,&argv);

	ostringstream os;
	os << "Code parallèle de calcul par différences finies" << endl << "Copyright EMSE, J. Bruchon" << endl;
	os << "Running on " << parallel::nb_proc() << " cores";
	parallel::print_message(os);

	double time1 = MPI_Wtime();
	
	// Initialisation step
	
	int dimension = 2;

	double point1[3] = {0,0,0};
	double point2[3] = {1,1,1};
	int np[3] = {300,300,30};
	
	grid une_grille(dimension,point1,point2,np);
	
	// Partition step
	
	partitioning my_partition(&une_grille,parallel::nb_proc());
 	my_partition.make_partition();
	une_grille.partitionne(my_partition.get_partition());
	my_partition.clean_partition();
	
	// Linear system definition and assembling step
	
	linear_system SL;
	int nb_col = 2*dimension + 1;
	SL.alloc_system(&une_grille,nb_col);
	
	double conductivite = 1;	
	assembling assembleur(&SL,&une_grille);
	bool CL_Dirichlet[4] = {1,1,1,1};
	double CL_D_values[4] = {1,2,3,4};
	assembleur.build_CL_Dirichlet(CL_Dirichlet,CL_D_values);
	assembleur.add_lapalacian(conductivite);	
	assembleur.apply_CL_Dirichlet();
	
	// Solver definition

	conjugate_gradient Solver;
	Solver.set_linear_system(&SL);
	Solver.solve_PCD();
	
	double time2 = MPI_Wtime();
	
	MPI_Barrier(parallel::comm);
	ostringstream os2;
	double dt = time2 - time1;
	os2 << "Temps total : " << dt;
	parallel::print_message(os2);
	
	// output step
	
	output* my_output=0;
	int* partition_field=0;
	double* global_solution=0;
	
	if ( parallel::IamMaster() )
	{
	  string output_name("partition");
	  my_output = new output_vtk();
	  my_output->set_filename(output_name);
	  
	  my_output->open_newfile();
	  my_output->write_grid(&une_grille);	
	}
	
	// have to reconstruct the field on master processor
	if ( !parallel::IamMaster() )
	{
	  int* local_part = une_grille.get_nodes();
	  MPI_Send(local_part, une_grille.get_NbNoeLocal(), MPI_INT, parallel::master, 0, parallel::comm);
	  double* sol = const_cast<double*>(SL.get_solution());
	  MPI_Send(sol, une_grille.get_NbNoeLocal(), MPI_DOUBLE, parallel::master, 1, parallel::comm);	  
	  delete [] local_part;
	}
	else
	{
	  int j;
	  partition_field = new int[une_grille.get_NbNoeGlobal()];
	  global_solution = new double[une_grille.get_NbNoeGlobal()];
	  
	  for (j=0; j<une_grille.get_NbNoeGlobal(); j++)
	    partition_field[j] = 0;
	 
	  for (int i=0; i<parallel::nb_proc(); i++)
	  {
	    if ( i != parallel::master )
	    {
	      int nbNoe_i = une_grille.get_NbNoePerProc(i);
	      int* local_part = new int[nbNoe_i];
	      MPI_Recv(local_part, nbNoe_i, MPI_INT, i, 0, parallel::comm, MPI_STATUS_IGNORE);
	      
	      for (j=0; j<nbNoe_i; j++)
	      {
		partition_field[local_part[j]] = i;		
	      }
	      
	      double* local_solution = new double[nbNoe_i];
	      MPI_Recv(local_solution, nbNoe_i, MPI_DOUBLE, i, 1, parallel::comm, MPI_STATUS_IGNORE);
	      
	     for (j=0; j<nbNoe_i; j++)
	      {
		global_solution[local_part[j]] = local_solution[j];		
	      }
	      
	      delete [] local_solution;
	      delete  [] local_part;	      
	    }
	  }
	  
	  int* local_part = une_grille.get_nodes();  // for master core
	  const double* px = SL.get_solution();
	  for (j=0; j<une_grille.get_NbNoeLocal(); j++)
	    global_solution[local_part[j]] = px[j];
	  delete [] local_part;
	  
	}
	
	
	if ( parallel::IamMaster() )
	{
	  my_output->write_header_point_data(une_grille.get_NbNoeGlobal());
	  my_output->write_scalarfield(une_grille.get_NbNoeGlobal(),"solution", global_solution);
	  my_output->write_scalarfield(une_grille.get_NbNoeGlobal(),"partition", partition_field);
	  my_output->close_file();
	}
	
	delete [] partition_field;
	delete [] global_solution;
	delete my_output;
		
	MPI_Finalize();
}



