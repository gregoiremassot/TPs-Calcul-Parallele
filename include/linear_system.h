#ifndef _LINEAR_SYSTEM_H_
#define _LINEAR_SYSTEM_H_

#include "grid.h"


/*
 * Classe permettant de gérer un système linéaire du type Ax = b avec un stockage Morse de la matrice A
 * (seules les entrées non nulles sont stockées) dans un environnement parallèle.
 */

class linear_system
{
public:
  linear_system() : A(0), b(0), x(0), une_grille(0) { ; }  // constructeur par défaut
  
  int alloc_system(const grid* lgrid, int lnb_col);
  virtual ~linear_system();   // destructeur pour désallouer la mémoire dans le tas 
  
  int matrix_vector_product(double* vecteur, double* res);
  
  
  const double* get_matrix() const { return A; }
  const double* get_right_member() const { return b; }
  const double* get_solution() const { return x; }
  
  const int& get_nb_lines() const { return nb_lines; }
  const int& get_nb_col() const { return nb_col; } 
  
  void set_matrix_to_zero();
  void set_right_member_to_zero();
  void set_solution_to_zero();
  
protected:
  double * A;
  double * b;
  double * x;
    
  const grid* une_grille;
  
  int nb_lines, nb_col;
};

#endif

