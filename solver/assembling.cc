
#include "assembling.h"
#include "grid.h"
#include "constants.h"
#include "parallel.h"

#include <math.h>
#include <iostream>
#include <mpi.h>

assembling::~assembling()
{
  delete [] CL_Dirichlet;
  delete [] CL_Neumann;
  delete [] CL_Dirichlet_Values;
  delete [] CL_Neumann_Values;
}

int assembling::build_CL_Dirichlet(const bool* CL, const double* values)
{
  int i,j;
  
  int nbNoe = grille->get_NbNoeLocal();
  int d = grille->get_space_dim();
  const double* corners = grille->get_corners();
  
  if ( CL_Dirichlet )
    delete CL_Dirichlet;
  CL_Dirichlet = new bool [nbNoe];
  
  if ( CL_Dirichlet_Values )
    delete CL_Dirichlet_Values;
  CL_Dirichlet_Values = new double [nbNoe];
  
  for (i=0; i<nbNoe; i++)
    CL_Dirichlet[i] = 0;
  
  for (i=0; i<nbNoe; i++)
  {
    const double* coord_i = grille->get_coordinates() + i*d;
    for (j=0; j<2*d; j++)
    {
      if ( CL[j] && ( fabs(coord_i[j%d] - corners[j]) < constants::precision ) )
      {
	CL_Dirichlet[i] = 1;
	CL_Dirichlet_Values[i] = values[j];
      }
    }
    
  }  
  
  return 0;
}

int assembling::apply_CL_Dirichlet()
{
  int i,j;
  int nbNoe = grille->get_NbNoeLocal();
  int nbcol = SL->get_nb_col();
  int D = grille->get_topological_dim();
  
  double* A = const_cast<double*>(SL->get_matrix());
  double* b = const_cast<double*>(SL->get_right_member());  
  
  if ( !CL_Dirichlet )
    return 1;
  
  // Etape 1 : envoi de l'information à chaque noeud
  int nb_proc = parallel::nb_proc();
  int rank = parallel::my_rank();
  vector<int>* collect_values_nodes;
  collect_values_nodes = new vector<int>[nb_proc];
  
  vector<double>* collect_values;
  collect_values = new vector<double>[nb_proc];
  
  for (i=0; i<nbNoe; i++)
  {
    if ( CL_Dirichlet[i] )
    {
      const int* pconnect = grille->get_connectivity(i);
      for (j=0; j<D; j++)
      {
	int voisin = pconnect[j];
	if ( voisin >= grille->get_NbNoeLocal() )
	{
	  int noe_glob = grille->get_node(i).first;
	  int proc = grille->get_node(voisin).second;
	  
	  collect_values_nodes[proc].push_back(noe_glob);
	  collect_values[proc].push_back(CL_Dirichlet_Values[i]);	  
	}	
      }   
    }  
  }
  
  // Envoi / reception de la taille des données
  for (i=0; i<nb_proc; i++)
  {
    int size = collect_values[i].size();
    if ( i != rank )
    {
      MPI_Send(&size, 1, MPI_INT, i, 3, parallel::comm);
    }
  }
  
  int* data_size = new int[nb_proc];
  for (i=0; i<nb_proc; i++)
  {
    if ( i != rank )
    {
      int* buff = data_size + i;
      MPI_Recv(buff, 1, MPI_INT, i, 3, parallel::comm,MPI_STATUS_IGNORE);     
    }
  }
  data_size[rank] = 0;
  
  // Envoi / réception des valeurs
  
  MPI_Request* request_send = new MPI_Request[2*nb_proc];
  MPI_Request* request_recv = new MPI_Request[2*nb_proc];
  MPI_Status* status = new MPI_Status[2*nb_proc];
  
  int count = 0;
  for (i=0; i<nb_proc; i++)
  {
    int size = collect_values[i].size();
    if ( ( i != rank ) && ( size > 0 ) )
    {
      MPI_Isend(collect_values_nodes[i].data(), size, MPI_INT, i, 0, parallel::comm,request_send+count++);      
      MPI_Isend(collect_values[i].data(), size, MPI_DOUBLE, i, 1, parallel::comm,request_send+count++);
    }
  }
  
  for (i=0; i<nb_proc; i++)
    collect_values[i].clear();
  
  for (i=0; i<nb_proc; i++)
    collect_values_nodes[i].clear();
   
  delete [] collect_values;
  delete [] collect_values_nodes;
  
  int total_size = 0;
  for (i=0; i<nb_proc; i++)
    total_size += data_size[i];
  
  int* node = new int[total_size];
  double* values = new double[total_size];
  
  int* p_node = node;
  double* p_values = values;   
  count = 0;
  for (i=0; i<nb_proc; i++)
  {
    int size = data_size[i];
    if ( ( i != rank ) && ( size > 0 ) )
    {
      MPI_Irecv(p_node, size, MPI_INT, i, 0, parallel::comm,request_recv+count++);
      MPI_Irecv(p_values, size, MPI_DOUBLE, i, 1, parallel::comm, request_recv+count++);
      
      p_node += size;
      p_values += size;
    }
  }
  
  MPI_Waitall(count,request_send,status);
  MPI_Waitall(count,request_recv,status);
  
  delete [] status;
  delete [] request_send;
  delete [] request_recv;
  
  map<int,double>* collect_pair = new map<int,double>[nb_proc];
  p_node = node;
  p_values = values; 
  for (i=0; i<nb_proc; i++)
  {
    int size = data_size[i];
    for (int u=0; u<size; u++)
    {
      collect_pair[i].insert(pair<int,double>(p_node[u],p_values[u])); 
    }
    p_node += size;
    p_values += size;
  }
  
  delete [] data_size;
  delete [] node;
  delete [] values;
  
  
  // Etape 2 : modification de la matrice et du second membre
  
  for (i=0; i<nbNoe; i++)
  {
    double* Ai = A + i*nbcol;
    
    if ( CL_Dirichlet[i] )
    {
      Ai[0] = 1.0;
      b[i] = CL_Dirichlet_Values[i];
      
      for (j=1; j<nbcol; j++)
	Ai[j] = 0;
    }
    else
    {
      const int* connect_i = grille->get_connectivity(i);
      for (j=1; j<nbcol; j++)
      {
	int node = connect_i[j-1];
	if ( ( node > -1 ) && ( node < nbNoe ) && CL_Dirichlet[node] )
	{	  
	    b[i] -= Ai[j]*CL_Dirichlet_Values[node];
	    Ai[j] = 0;	  
	}
	else if ( node >= nbNoe )
	{
	  int num_glob = grille->get_node(node).first;
	  int proc = grille->get_node(node).second;
	  
	  map<int,double>::const_iterator iter = collect_pair[proc].find(num_glob);
	  if ( iter != collect_pair[proc].end() )
	  { 
	    b[i] -= Ai[j]*iter->second;
	    Ai[j] = 0;
	  }	  	  
	}
      }      
    }
  }
  
  return 0;
}



int assembling::add_lapalacian(double k)
{
  int nbNoe = grille->get_NbNoeLocal();
  int d = grille->get_space_dim();
  int nbcol = SL->get_nb_col();
  
  double* A = const_cast<double*>(SL->get_matrix());
  double* b = const_cast<double*>(SL->get_right_member());
  
  int i,j;
  for (i=0; i<nbNoe; i++)
  {
    double* Ai = A + i*nbcol;
    for (j=0; j<d; j++)
    {
      double pas_spatial = grille->get_dxdydz(j);
      double lambda = k/(pas_spatial*pas_spatial);
      
      Ai[0] += 2.0*lambda;
      Ai[2*j+1] += -lambda;
      Ai[2*j+2] += -lambda;           
    }
    
 
   /* std::cout << "ligne " << i << " : ";
    for (j=0; j<nbcol; j++)
      std::cout << Ai[j] << "  ";
    std::cout << std::endl;      */
  }
  

  return 0;
}
