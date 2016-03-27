#ifndef _PARTITIONING_H_
#define _PARTITIONING_H_

#include "grid.h"
#include <vector>

class partitioning
{
public:
	partitioning() : nb_domains(0), partition(0), agrid(0) { ; }
	partitioning(grid* lgrid, int lnb_domains);
	~partitioning();

	int make_partition();
	int clean_partition();

	void set_grid(grid* lgrid) { agrid = lgrid; d = agrid->get_space_dim(); }
	void set_nbdomains(int lnb_domains) {nb_domains = lnb_domains; }

	const int& get_nbdomains() const { return nb_domains; }
	const int* get_partition() const { return partition; }

protected:
	int nb_domains;
	int* partition;
	grid* agrid;

	int d;

	// Functions for building the partition of a mesh (maybe do a derivated class?
	int grid_partition();
	int first_node(const int* degree);
	int find_nucleus(const vector<int>& front, const int* node_state, int nbNoe);
};


#endif
