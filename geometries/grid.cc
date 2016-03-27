
#include <iostream>
#include <sstream>
#include <math.h>
#include <set>

#include "grid.h"
#include "parallel.h"

grid::grid()
{
	initialise();
}

grid::~grid()
{
  delete [] corners;
  delete [] connectivity;
  delete [] coordinates;
  delete [] nbNoe_per_proc;
  nodes.clear();
}

void grid::initialise()
{
	D = 0; d = 0;
	nbNoe = 0; nbNoe_loc = 0; nbNoe_glob = 0;

	connectivity = 0; coordinates = 0;
	nbNoe_per_proc = 0;
	corners = 0;
}

grid::grid(string & grid_file)
{
	 initialise();
	 ifstream afile(grid_file.c_str(), ios::in);
	 int err = 0;

	 if (afile)
	 {
		err = read_file(afile);
		afile.close();

		if ( err )
		{
			ostringstream os;
			os << "The format of file " << grid_file << " is not correct!";
			parallel::print_problem(os);
		}
	 }
	 else
	 {
		 ostringstream os;
		 os << "File " << grid_file << " unknown!";
		 parallel::print_problem(os);
	 }
}

int* grid::get_nodes(int* my_nodes)
{
  if ( !my_nodes )
    my_nodes = new int[nbNoe_loc];
  
  int i;
  vector<pair<int,int> >::const_iterator it;
  for (it=nodes.begin(), i=0; i<nbNoe_loc; i++, it++)
  {
    my_nodes[i] = it->first;
  }
  
  return my_nodes;
}

int grid::read_file(ifstream& afile)
{
	int i;
	int err  = 0;

	parallel::print_message("Grid reading in progress...");
	afile >> d;
	D = 2*d;
	nbNoe = 1;
	
	for (i=0; i<d; i++)
	{
	  afile >> nbNoe_xyz[i];
	  nbNoe *= nbNoe_xyz[i];
	  if ( nbNoe <= 0 )
	    err = 1;
	}

	for (i=d; i<3; i++)
	  nbNoe_xyz[i] = 1;
	
	for (i=0; i<d; i++)
	{
	  afile >> dxdydz[i];
	  if ( dxdydz[i] <= 0 )
	    err = 1;
	}

	for (i=d; i<3; i++)
	  dxdydz[i] = 1;

	if ( ( d <= 0 ) || err )
	  return 1;

	coordinates = new double [d*nbNoe];
	connectivity = new int [D*nbNoe];

	for (i=0; i<d*nbNoe; i++)
	  afile >> coordinates[i];

	for (i=0; i<D*nbNoe; i++)
	  afile >> connectivity[i];

	nbNoe_loc = nbNoe;
	nbNoe_glob = nbNoe;
	nbNoe_per_proc = new int[1];
	nbNoe_per_proc[0] = nbNoe;
	
	ostringstream lstr;
	lstr << "Space dimension: " << d << " ; topological dimension: " << D << endl;
	lstr << "Number of nodes: " << nbNoe << endl;
	lstr << "Reading step done";
	
	// Missing : calculation of corners
	parallel::print_problem("Dans grid.cc : calcul du tableau corners non encore implémenté.");
	return err;
}

grid::grid(int dim, double* corner1, double* corner2, int* nbxyz)
{
  int i;
  int err = 0;
  parallel::print_message("Grid construction in progress...");
  initialise();
  d = dim;
  D = 2*d;
    
  nbNoe = 1;
  for (i=0; i<d; i++)
  {
    nbNoe_xyz[i] = nbxyz[i];
    nbNoe *= nbNoe_xyz[i];
    if ( nbNoe <= 0 )
      err = 1;     
  }

  for (i=d; i<3; i++)
    nbNoe_xyz[i] = 1;
	  
  coordinates = new double [d*nbNoe];
  connectivity = new int [D*nbNoe];
  
  for (i=0; i<d; i++)
  {
    dxdydz[i] = fabs(corner1[i] - corner2[i])/((double) nbxyz[i]-1.0);
    if ( dxdydz[i] <= 0 )
      err = 1;
  }
  
  for (i=d; i<3; i++)
    dxdydz[i] = 1;
  
  if ( ( d <= 0 ) || err )
    parallel::print_problem("Bad arguments in the grid constructor");
  
  double corner_inf[3];
  double corner_sup[3];
  double coord[3];
  
  for (i=0; i<d; i++)
  {
    if ( corner1[i] < corner2[i] )
    {
      corner_inf[i] = corner1[i];
      corner_sup[i] = corner2[i];
    }
    else
    {
      corner_inf[i] = corner2[i];
      corner_sup[i] = corner1[i];
    }
  }
  
  for (i=d; i<3; i++)
  {
    corner_inf[i] = 0;
    corner_sup[i] = 0;
  }
    
  for (i=0; i<3; i++)
    coord[i] = corner_inf[i];
  
  int count = 0;
  for (int nz=0; nz<nbNoe_xyz[2]; nz++)
  {
    for (int ny=0; ny<nbNoe_xyz[1]; ny++)
    {
      for (int nx=0; nx<nbNoe_xyz[0]; nx++)
      {
	for (i=0; i<d; i++)
	  coordinates[count*d+i] = coord[i];

	count++;
	coord[0] += dxdydz[0];
	if ( ( coord[0] > corner_sup[0] ) || ( corner_sup[0] - coord[0] < 0.5*dxdydz[0] ) )
	  coord[0] = corner_sup[0];
      }
      coord[0] = corner_inf[0];
      coord[1] += dxdydz[1];
      if ( ( coord[1] > corner_sup[1] ) || ( corner_sup[1] - coord[1] < 0.5*dxdydz[1] ) )
	  coord[1] = corner_sup[1];
    }
    coord[1] = corner_inf[1];
    coord[2] += dxdydz[2];
    if ( ( coord[2] > corner_sup[2] ) || ( corner_sup[2] - coord[2] < 0.5*dxdydz[2] ) )
	  coord[2] = corner_sup[2];
  }
  
  // Construction of the connectivity : x-1, x+1, y-1, y+1, z-1, z+1
  
  for (i=0; i<nbNoe; i++)
  {
    int shift = 1;
    int* p_connec = connectivity + i*D;
    const double* p_coord = get_coordinates(i);
    for (int k=0; k<d; k++)
    {
      int left = i - shift;
      int right = i + shift;

      if ( fabs(p_coord[k] - corner_inf[k]) < 0.1*dxdydz[k] )
	left = -1;
      
      if ( fabs(p_coord[k] - corner_sup[k]) < 0.1*dxdydz[k] )
	right = -1;
	
      p_connec[2*k] = left;
      p_connec[2*k+1] = right;
      shift *= nbNoe_xyz[k];
    }
  }
      
  nbNoe_loc = nbNoe;
  nbNoe_glob = nbNoe;
  nbNoe_per_proc = new int[1];
  nbNoe_per_proc[0] = nbNoe;
  
  corners = new double[2*d];
  for (i=0; i<d; i++)
    corners[i] = corner_inf[i];
  
  for (i=0; i<d; i++)
    corners[i+d] = corner_sup[i];
  
  ostringstream lstr;
  lstr << "Space dimension: " << d << " ; topological dimension: " << D << endl;
  lstr << "Number of nodes: " << nbNoe << endl;
  lstr << "Grid construction done";  
  parallel::print_message(lstr);
}

int grid::partitionne(const int* partition)
{
  int i,j;
  const int* p_connec;
  
  int moi = parallel::my_rank();
  nbNoe_loc = 0;
 
  // Step 1. Searching the nodes for processor "moi", including the nodes of closest neighbourhood
  set<int> local_nodes;
  for (i=0; i<nbNoe; i++)
  {
    if ( partition[i] == moi )
    {
      nbNoe_loc++;
      local_nodes.insert(i);
      p_connec = get_connectivity(i);
      for (j=0; j<D; j++)
      {
	int neigh = p_connec[j];
	if ( ( neigh != -1 ) && ( partition[neigh] != moi ) )
	  local_nodes.insert(neigh);
      }
    }
  }
  
  nbNoe = local_nodes.size();  // number of nodes that belong to processor "moi"
  
  
  // Stepz 2. Construction of the "nodes" container
  set<int>::const_iterator it_set;
  nodes.clear();
  for (it_set = local_nodes.begin(); it_set != local_nodes.end(); it_set++)
  {
    if ( partition[*it_set] == moi )
      nodes.push_back(pair<int,int>(*it_set,moi));
  }
  
  for (it_set = local_nodes.begin(); it_set != local_nodes.end(); it_set++)
  {
    if ( partition[*it_set] != moi )
      nodes.push_back(pair<int,int>(*it_set,partition[*it_set]));
  }
  local_nodes.clear();
  
  // Step 3. Construction of the local array of coordinnates
  double * local_coordinates = new double[nbNoe_loc*d];
  vector<pair<int,int> >::const_iterator it_list;
  for (it_list=nodes.begin(), i=0; i<nbNoe_loc; it_list++, i++)
  {
    int noe_glob = it_list->first;
    double* p_coord_loc = local_coordinates + i*d;
    for (j=0; j<d; j++)
    {
      p_coord_loc[j] = get_coordinates(noe_glob)[j];
    }    
  } 
  
  delete [] coordinates;
  coordinates = local_coordinates;
  
  // Step 4. Construction of the local array of connectivities
  int* local_connectivity = new int[nbNoe_loc*D];
  int* glob_to_loc = new int[nbNoe_glob];
  for (it_list=nodes.begin(), i=0; it_list!=nodes.end(); it_list++, i++)
  {
    int noe_glob = it_list->first;
    glob_to_loc[noe_glob] = i;
  }  
  
  for (it_list=nodes.begin(), i=0; i<nbNoe_loc; it_list++, i++)
  {
    int noe_glob = it_list->first;
    int* p_connec_loc = local_connectivity + i*D;
    for (j=0; j<D; j++)
    {
      p_connec_loc[j] = glob_to_loc[get_connectivity(noe_glob)[j]];
    }    
  } 
  
  delete [] glob_to_loc;
  delete [] connectivity;
  connectivity = local_connectivity;
  
  // Step 5. Each proc know how many nodes are treated by each processor
  delete [] nbNoe_per_proc;
  nbNoe_per_proc = new int[parallel::nb_proc()];
  MPI_Allgather(&nbNoe_loc,1,MPI_INT,nbNoe_per_proc,1,MPI_INT,parallel::comm);
    
  return 0;
}