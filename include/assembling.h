#ifndef _ASSEMBLING_H_
#define _ASSEMBLING_H_

#include "linear_system.h"


/*
 * Classe permettant d'assembler un système linéaire issu d'un discrétisation par DF
 */

class assembling
{
public:
  assembling(linear_system* lSL=0, const grid* lgrid=0) : SL(lSL), grille(lgrid), CL_Dirichlet(0), CL_Neumann(0),
							  CL_Dirichlet_Values(0), CL_Neumann_Values(0) { ; }  // constructeur par défaut
  
  virtual ~assembling();  // destructeur pour désallouer la mémoire dans le tas 
  void set_linear_system(linear_system* lSL) { SL = lSL; }
  void set_grid(const grid* lgrille) { grille = lgrille; }
  
  int build_CL_Dirichlet(const bool* CL, const double* values);
  int build_CL_Neumannt(const bool* CL, const double* values) { ; }
  
  int apply_CL_Dirichlet();
  
  int add_lapalacian(double k);
  
protected:
  linear_system* SL;
  const grid* grille;
  
  bool* CL_Dirichlet;
  bool* CL_Neumann;
  double* CL_Dirichlet_Values;
  double* CL_Neumann_Values;
};

#endif

