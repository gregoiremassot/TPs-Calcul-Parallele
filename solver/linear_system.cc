#include "linear_system.h"
#include "parallel.h"
#include <list>
#include <map>

using namespace std;

linear_system::~linear_system()
{
  delete[] A;
  delete[] b;
  delete[] x;
}


int linear_system::alloc_system(const grid* lgrid, int lnb_col)
{
  une_grille = lgrid;
  nb_lines = une_grille->get_NbNoeLocal();
  nb_col = lnb_col;
  
  if ( A )
    delete [] A;
  
  A = new double [nb_lines*nb_col];

  if ( b )
    delete [] b;
  
  b = new double [nb_lines];
  

  if ( x )
    delete [] x;
  
  x = new double [nb_lines];
  
  set_matrix_to_zero();
  set_right_member_to_zero();
  set_solution_to_zero();
  
  return 0;
}

void linear_system::set_matrix_to_zero()
{
 for (int i=0; i<nb_lines*nb_col; i++)
   A[i] = 0;   
}

void linear_system::set_right_member_to_zero()
{
 for (int i=0; i<nb_lines; i++)
   b[i] = 0;   
}

void linear_system::set_solution_to_zero()
{
 for (int i=0; i<nb_lines; i++)
   x[i] = 0;   
}

int linear_system::matrix_vector_product(double* vecteur, double* res)
{
  int i,j;
  int nb_proc = parallel::nb_proc();
  int rank = parallel::my_rank();
  int D = une_grille->get_topological_dim();
  
  if ( !res )
    res = new double[nb_lines];
  
  // Step 1.  Chaque coeur envoie aux coeurs voisins les valeurs de vector dont ils ont besoin, et reçoit les siennes
  
  vector<int>* collect_values_nodes;
  collect_values_nodes = new vector<int>[nb_proc];
  
  vector<double>* collect_values;
  collect_values = new vector<double>[nb_proc];
  
  for (i=0; i<une_grille->get_NbNoeLocal(); i++)
  {
    const int* pconnect = une_grille->get_connectivity(i);
    for (j=0; j<D; j++)
    {
      int voisin = pconnect[j];
      if ( voisin >= une_grille->get_NbNoeLocal() )
      {
	int noe_glob = une_grille->get_node(i).first;
	int proc = une_grille->get_node(voisin).second;
	
	collect_values_nodes[proc].push_back(noe_glob);
	collect_values[proc].push_back(vecteur[i]);
      }
    }   
  }
    int total_size = 0;
  for (i=0; i<nb_proc; i++)
    total_size += collect_values[i].size();
  
  int* data_size = new int[nb_proc];
  for (i=0; i<nb_proc; i++)
    data_size[i] = collect_values[i].size();
    
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
    collect_values_nodes[i].clear();
  delete [] collect_values_nodes;
  
  for (i=0; i<nb_proc; i++)
    collect_values[i].clear();
  delete [] collect_values;

  int* node = new int[total_size];
  int* p_node = node;
  double* values = new double[total_size];
  double* p_values = values;  
   
    // Cas particulier des grilles structurées : un coeur envoie et reçoit la même quantité d'infos. Les vecteurs sont donc déjà bien dimensionnés
  // On écrase donc les données envoyées avec celles reçues
  
  count  = 0;
  for (i=0; i<nb_proc; i++)
  {
    int size = data_size[i];
    if ( ( i != rank ) && ( size > 0 ) )
    {
      MPI_Irecv(p_node, size, MPI_INT, i, 0, parallel::comm,request_recv+count++);
      p_node += size;
      
      MPI_Irecv(p_values, size, MPI_DOUBLE, i, 1, parallel::comm, request_recv+count++);
      p_values += size;
    }
  }
  
  
  MPI_Waitall(count,request_send,status);
  MPI_Waitall(count,request_recv,status);
    
  delete [] status;
  delete [] request_send;
  delete [] request_recv;  
  
  p_node = node;
  p_values = values;
  
  map<int,double>* collect_pair = new map<int,double>[nb_proc];
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
  
  // Step 2. Matrix - vector product 
  
  // Le tableau des connectivités sert ici à connaître à quel noeud de la grille locale correspond la colonne j, avec 0 <= j < nb_col
  // On suppose que le nombre de valeurs non nulle sur une ligne i n'excède pas 1 + nombre de voisins du noeud i.
  // On suppose que l'élément diagonal de la ligne i est toujours stocké en première position.
  
  for (i=0; i<nb_lines; i++)
  {
    double* p_Ai = A + i*nb_col;  // ligne i de la matrice
    res[i] = p_Ai[0]*vecteur[i];
    
    const int* connect_pi = une_grille->get_connectivity(i);
    for (j=0; j<nb_col-1; j++)
    {
      int colonne_j = connect_pi[j];
      double val = 0;
      if ( ( colonne_j > -1 ) && ( colonne_j < une_grille->get_NbNoeLocal() ) )
      {
	val =  vecteur[colonne_j];
      }
      else if ( colonne_j > -1 )
      {
	int noej_glob = une_grille->get_node(colonne_j).first;
	int proc = une_grille->get_node(colonne_j).second;
	
	map<int,double>::const_iterator iter = collect_pair[proc].find(noej_glob);
	if ( iter == collect_pair[proc].end() )
	  parallel::print_problem("Problème dans le produit matrice - vecteur de linear_system.cc");
	else
	  val = iter->second;
      }
      
      res[i] += p_Ai[j+1]*val;
    }    
  }
  
   
  for (i=0; i<nb_proc; i++)
    collect_pair[i].clear();
   
  delete [] collect_pair;
  
  return 0;
}
