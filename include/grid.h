#ifndef _GRID_H_
#define _GRID_H_

#include <fstream>
#include <map>
#include <vector>

using namespace std;

/*
 * Classe permettant de gérer une grille cartésienne. Les données membres de la classe sont :
 * 	- le tableau de coordonnées des noeuds, double* coordinates;
 * 	- le tableau des coonectivités des noeuds, int* node_connectivity;
 * 	- le tableau contenant le numéro global des noeuds (calcul parallèle) et le processeur auquel ils appartiennent (noeuds frontière)
 * 	  map<int,int>* nodes;
 * 
 * On accède à ces données depuis l'extérieur par des fonctions d'accès get_...(...)
 * 
 */

class grid
{
public:
	grid(); 		  // constructeur par défaut
	grid(string & grid_file); // constructeur à partir du fichier contenant la grille
	grid(int dim, double* corner1, double* corner2, int* nbxyz);   // constructeur pour une grille "parallélépipédique" 
	int partitionne(const int* partition);
	
	virtual ~grid();          // destructeur pour désallouer la mémoire dans le tas 

	const int& get_space_dim() const { return d; }
	const int& get_topological_dim() const { return D; }
	const int& get_nbNoe() const { return nbNoe; }
	const int& get_NbNoeGlobal() const { return nbNoe_glob; }
	const int& get_NbNoeLocal() const { return nbNoe_loc; }
	const int& get_NbNoePerProc(int i)  { return nbNoe_per_proc[i]; }


	const double* get_corners() const { return corners; }
	const double* get_coordinates() const { return coordinates; }
	const double* get_coordinates(int i) const { return coordinates + i*d; }
	
	const int* get_connectivity() const { return connectivity; }
	const int* get_connectivity(int i) const { return connectivity + i*D; }
	
	int* get_nodes(int* my_nodes=0);
	const pair<int,int>& get_node(int i) const { return nodes[i]; }
	
	const double& get_nbNoe_xyz(int i) const { return nbNoe_xyz[i]; }
	const double& get_dxdydz(int i) const { return dxdydz[i]; }

protected:

	void initialise();
	int read_file(ifstream& afile);

// Scalar parameters
	int d,D;   	// d : dimension spatiale, D : nombre de voisins (dimension topologique)
	int nbNoe;	// nombre de noeuds total du domaine de la grille sur le processeur
	int nbNoe_loc;  // nombre de noeuds sur le processeur local, sans compter les noeuds voisins du processeur
	int nbNoe_glob; // nombre de noeuds total
	

// Arrays
	double* corners;
	double * coordinates;
	int * connectivity;
	vector< pair<int,int> > nodes;
	
	double nbNoe_xyz[3];  // nombre de noeuds selon x, y, z, sert pour l'affichage
	double dxdydz[3];     // pas d'espace selon x, y et z, ici constrant pour simplifier l'affichage	
	
// For parallel computing
	int* nbNoe_per_proc;
};

#endif
