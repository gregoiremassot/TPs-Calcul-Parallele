
#include <iostream>
#include <stdlib.h>
#include <set>
#include <vector>
#include <algorithm>
#include <mpi.h>
#include "math.h"

#include "partitioning.h"
#include "parallel.h"
#include "constants.h"


using namespace std;

template<typename T1, typename T2> 
bool pair_comparison(const pair<T1, T2> & i, const pair<T1, T2> & j) { return (i.second < j.second); }

partitioning::partitioning(grid* lgrid, int lnb_domains)
{
	set_grid(lgrid);
	set_nbdomains(lnb_domains);
	partition = 0;
}

partitioning::~partitioning()
{
  delete [] partition;
}

int partitioning::make_partition()
{
	if ( !agrid )
		return 1;

	return grid_partition();
}

int partitioning::clean_partition()
{
  if ( partition )
    delete [] partition;
  partition = 0;
}

int partitioning::grid_partition()
{
	int D = agrid->get_topological_dim();
	int nbNoe = agrid->get_nbNoe();
	int i,j;
	const int* p_connect;
	
	// Dynamic allocation of requiered data arrays
	
	if ( partition )
	  delete [] partition;
	
	partition = new int[nbNoe];
	int* current_degree = new int[nbNoe];
	int* local_nbNoe = new int[nb_domains];

	// Initialisation of these arrays
	for (i=0; i<nbNoe; i++)
	  partition[i] = 0;

	for (i=0; i<nbNoe; i++)
	{
	 current_degree[i] = D;
	 p_connect =  agrid->get_connectivity(i);
	 for (j=0; j<D; j++)
	 {
	   if ( p_connect[j] == -1 )
	     current_degree[i]--;
	 }
	}
	
	int nucleus = first_node(current_degree);
	vector<int> generation_1;
	set<int> generation_2;
	
	int sum_ni = 0;

	vector<int> old_front;
	for (i=1; i<nb_domains; i++)
	{
	  int cur_nbNoe = 0;
	  generation_1.push_back(nucleus);
	  int nbNoe_loc = (nbNoe - sum_ni)/(nb_domains - i+1);
	  	  
	  bool again = 1;
	  while ( ( cur_nbNoe < nbNoe_loc ) && again )
	  {
	    int z;
	    again = 0;
	    
	    // look at if we can take all the unmarked nodes in generation 1, or if we have to choose those having a minimal current degree
	    
	    vector<int>::const_iterator it;
	    int K = (int) generation_1.size();
	    
	    if ( cur_nbNoe + K > nbNoe_loc )   // have to take only the nbNoe_loc - cur_nbNoe first nodes sorted by minimal current degree
	    {
	      vector< pair<int,int> > generation1_bis(generation_1.size());
	      for (z=0, it=generation_1.begin(); it!=generation_1.end(); z++, it++)
	      {
		generation1_bis[z] = pair<int,int>(*it,current_degree[*it]);
	      }
	      
	      sort(generation1_bis.begin(),generation1_bis.end(),pair_comparison<int,int>);  // sort by minimal value of current degree
	      	      
	      vector< pair<int,int> >::const_iterator it_vec = generation1_bis.begin();
	      for (z=0; z<generation_1.size(); z++, it_vec++)
	      {
		generation_1[z] = it_vec->first;		
	      }	      
	    }
	    
	    it = generation_1.begin();
	    while ( ( it != generation_1.end() ) && ( cur_nbNoe < nbNoe_loc ) )
	    {
	      partition[*it] = i;
	      cur_nbNoe++;
	      again = 1;
	      
	      p_connect = agrid->get_connectivity(*it);
	      for (j=0; j<D; j++)
	      {
		int neighbour = p_connect[j];

		if ( ( neighbour != -1 ) && !partition[neighbour] )
		{
		  generation_2.insert(neighbour);
		  current_degree[neighbour]--;		  
		}		
	      }
	      it++;	      
	    }
	    
	    while ( it != generation_1.end() )
	    {
	      generation_2.insert(*it);
	      it++;
	    }
	    
	    generation_1.clear();
	    for (set<int>::const_iterator its = generation_2.begin(); its != generation_2.end(); its++)
	    {
		generation_1.push_back(*its);
	    }
	    
	    generation_2.clear();
	    if ( !again  && (fabs(nbNoe_loc - cur_nbNoe) > 0.4*nbNoe_loc) )
	    {    // means that we are in a local dead end, so have to remove our last actions
	      again = 1;
	      int old_dom = 0;
	      
	      int z;
	      int * dom = new int[nb_domains+1];
	      for (z=0; z<=nb_domains; z++)
		dom[z] = 0;
	      
	      p_connect = agrid->get_connectivity(nucleus);
	      for (j=0; j<D; j++)
	      {
		int neighbour = p_connect[j];
		if ( neighbour != -1 )
		{
		  int u = partition[neighbour];
		  dom[u]++;
		}
	      }
	      
	      int max = 0;
	      for (z=1; z<=nb_domains; z++)
	      {
		if ( ( z != i ) && ( dom[z] > max ) )
		{
		  max = dom[z];
		  old_dom = z;
		}
	      }
	      delete [] dom;
	      {
		ostringstream os;
		os << cur_nbNoe << " nodes are reallocated to domain number " << old_dom;
		parallel::print_message(os);
	      }
	      
	      for (z=0; z<nbNoe; z++)
	      {
		if ( partition[z] == i )
		{
		  partition[z] = old_dom;	
		}		
	      }
	      
	      if ( old_dom )
		local_nbNoe[old_dom-1] += cur_nbNoe;
	      
	      sum_ni += cur_nbNoe;
	      nbNoe_loc = (nbNoe - sum_ni)/(nb_domains - i+1);
	      cur_nbNoe = 0;
	      
	      generation_1.clear();
	      nucleus = find_nucleus(old_front,partition,nbNoe);
	      for (vector<int>::iterator itv = old_front.begin(); itv != old_front.end(); itv++)
		if ( *itv == nucleus )
		  old_front.erase(itv);
		
		generation_1.push_back(nucleus);	    
	    }
	  }
	  
	  local_nbNoe[i] = cur_nbNoe;
	  sum_ni += cur_nbNoe;
	  
	  // Search next nucleus
	  
	  if ( i != nb_domains - 1 )
	  {
	    old_front.clear();
	    
	    vector< pair<int,int> > generation1_bis;
	    vector<int>::const_iterator	it = generation_1.begin();
	    while ( it != generation_1.end() )
	    {
	      generation1_bis.push_back(pair<int,int>(*it,current_degree[*it]));
	      it++;	      
	    }
	    sort(generation1_bis.begin(),generation1_bis.end(),pair_comparison<int,int>);  // sort by minimal value of current degree
	    vector< pair<int,int> >::const_iterator it_vec = generation1_bis.begin();
	    while ( it_vec != generation1_bis.end() )
	    {
	      old_front.push_back(it_vec->first);
	      it_vec++;
	    }
	    generation_1.clear();
	    nucleus = find_nucleus(old_front,partition,nbNoe);
	    for (vector<int>::iterator itv = old_front.begin(); itv != old_front.end(); itv++)
	      if ( *itv == nucleus )
		old_front.erase(itv);	  
	  }	
	}

	
	local_nbNoe[0] = nbNoe - sum_ni;

	for (i=0; i<nb_domains; i++)
	{
	  ostringstream os;
	  os << "Domain " << i << " contains " << local_nbNoe[i] << " nodes";
	  parallel::print_message(os);
	}

	generation_1.clear();
	generation_2.clear();

	delete [] local_nbNoe;
	delete [] current_degree;

	return 0;
}

int partitioning::find_nucleus(const vector<int>& front, const int* node_state, int nbNoe)
{
	int nucleus = 0;
	vector<int>::const_iterator it = front.begin();
	while ( ( it != front.end() ) && !nucleus )
	{
		if ( !node_state[*it] )
			nucleus = *it;
		it++;
	}
	if ( !nucleus )
	{
		int z = 0;
		while ( !nucleus && ( z <= nbNoe ) )
		{
			if ( !node_state[z] )
				nucleus = z;
			z++;
		}
		if ( !nucleus )
		{
			ostringstream os;
			os << "In partitioning class: can not split the grid in " << nb_domains << " subdomains!";
			parallel::print_problem(os);
		}
	}
	return nucleus;
}

int partitioning::first_node(const int* degree)
{
  int node = 0;
  int i;
  int nbNoe = agrid->get_nbNoe();
  int connect = agrid->get_topological_dim();
 
  for (i=0; i<nbNoe; i++)
  {
    if ( degree[i] < connect )
    {
      node = i;
      connect = degree[i];
    }
  }
  return node;
}
