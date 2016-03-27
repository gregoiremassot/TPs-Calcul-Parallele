#include "conjugate_gradient.h"
#include "constants.h"
#include "parallel.h"

using namespace std;

conjugate_gradient::conjugate_gradient()
{
  SL = 0;
  precision = 1.0e-6; // constants::precision;
  nbItermax = 10000;  
}

double conjugate_gradient::inner_product(const double* vec1, const double* vec2) const
{
  int nb_lines = SL->get_nb_lines();
  double ps = 0;
  for (int i=0; i<nb_lines; i++)
    ps += vec1[i]*vec2[i];
  
  double somme_totale = 0;
  
  MPI_Allreduce(&ps, &somme_totale, 1, MPI_DOUBLE, MPI_SUM, parallel::comm);

  return somme_totale;
}

int conjugate_gradient::solve()
{
  int j;
  
  int nb_lines = SL->get_nb_lines();
  double* r = new double[nb_lines];   // gradient de -J = -((Ax,x)/2 - (b,x)) = residu b - Ax
  double* descente = new double[nb_lines];
  double* A_fois_descente = new double[nb_lines];
  
  double* x = const_cast<double*>(SL->get_solution());
  const double* b = SL->get_right_member();
  
  double r_i_fois_r_i = 0;
 
  // Initialisation
  
  SL->matrix_vector_product(x,r);
  for (j=0; j<nb_lines; j++)
  {
    r[j] = b[j] - r[j];
  }
  for (j=0; j<nb_lines; j++)
    descente[j] = r[j];
  
  r_i_fois_r_i = inner_product(r,r);
  
  int i = 0;
  bool go = 1;
  
  if ( r_i_fois_r_i < precision )
    go = 0;
  
  ostringstream os;
  os << "GC - Résidu initial : " <<  std::scientific << r_i_fois_r_i << endl;
  parallel::print_message(os);	
  
  // Boucle principale
	
  while ( go )
  {
    double alpha_i = inner_product(r,descente);
    SL->matrix_vector_product(descente,A_fois_descente);
    alpha_i /= inner_product(descente,A_fois_descente);
    
    for (j=0; j<nb_lines; j++)
      x[j] += alpha_i*descente[j];
    
    for (j=0; j<nb_lines; j++)
      r[j] -= alpha_i*A_fois_descente[j];
    
    double r_ip1_2 = inner_product(r,r);
    
    if ( r_ip1_2 > precision )
    {
      double beta_i = r_ip1_2/r_i_fois_r_i;
      for (j=0; j<nb_lines; j++)
	descente[j] = r[j] + beta_i*descente[j];
      
      r_i_fois_r_i = r_ip1_2;      
    }
    else
      go = 0;
    
    i++;
    
    if ( i > nbItermax )
      go = 0;
    
    if  ( !(i %50) || !go )
    {
        ostringstream os2;
	os2 << "GC - itértation " << i << ", résidu = " <<  std::scientific << r_ip1_2 << endl;
	parallel::print_message(os2);
    }
    
  }
    
  delete [] r;
  delete [] descente;
  delete [] A_fois_descente;
}


int conjugate_gradient::solve_PCD()
{
  int j;
  
  int nb_lines = SL->get_nb_lines();
  double* r = new double[nb_lines];   // gradient de -J = -((Ax,x)/2 - (b,x)) = residu b - Ax
  double* z = new double[nb_lines];
  double* descente = new double[nb_lines];
  double* A_fois_descente = new double[nb_lines];
  
  double* x = const_cast<double*>(SL->get_solution());
  const double* b = SL->get_right_member();
  const double* A = SL->get_matrix();
  
  double r_i_fois_r_i = 0;
  int nbcol = SL->get_nb_col();
 
  // Initialisation
  SL->matrix_vector_product(x,r);
  for (j=0; j<nb_lines; j++)
  {
    r[j] = b[j] - r[j];
  }
   
  for (j=0; j<nb_lines; j++)
  {
    z[j] = r[j]/A[j*nbcol];
  }
  
  for (j=0; j<nb_lines; j++)
    descente[j] = z[j];
  
  r_i_fois_r_i = inner_product(r,r);
  double zi_ri = inner_product(z,r);
  
  int i = 0;
  bool go = 1;
  
  if ( r_i_fois_r_i < precision )
    go = 0;
  
  ostringstream os;
  os << "GCP Diag - Résidu initial : " <<  std::scientific << r_i_fois_r_i;
  parallel::print_message(os);	
    
  // Boucle principale
	
  while ( go )
  {
    double alpha_i = zi_ri;
    SL->matrix_vector_product(descente,A_fois_descente);
    alpha_i /= inner_product(descente,A_fois_descente);
    
    for (j=0; j<nb_lines; j++)
      x[j] += alpha_i*descente[j];
    
    for (j=0; j<nb_lines; j++)
      r[j] -= alpha_i*A_fois_descente[j];
    
    double r_ip1_2 = inner_product(r,r);
    
    if ( r_ip1_2 > precision )
    {
      for (j=0; j<nb_lines; j++)
      {    
	z[j] = r[j]/A[j*nbcol];	
      }
      double zi_ri_new = inner_product(r,z);
      
      double beta_i = zi_ri_new/zi_ri;
      for (j=0; j<nb_lines; j++)
	descente[j] = z[j] + beta_i*descente[j];
      
      r_i_fois_r_i = r_ip1_2;     
      zi_ri = zi_ri_new;
    }
    else
      go = 0;
    
    i++;
    
    if ( i > nbItermax )
      go = 0;
    
    if  ( !(i %50) || !go )
    {
        ostringstream os2;
	os2 << "GCP Diag - itértation " << i << ", résidu = " <<  std::scientific << r_ip1_2;
	parallel::print_message(os2);
    }
    
  }
    
  delete [] r;
  delete [] z;
  delete [] descente;
  delete [] A_fois_descente;
}