
#include "output_vtk.h"
#include "grid.h"
#include "parallel.h"

#include <iostream>

const string output_vtk::extension = "vtk";

const string& output_vtk::get_extension() const
{
	return output_vtk::extension;
}

void output_vtk::write_header_point_data(int nbV)
{
  *my_file << "POINT_DATA "<< nbV << endl;  // write this only the first time
}

void output_vtk::write_header()
{
	*my_file << "# vtk DataFile Version 2.0" << endl;
	*my_file << full_name << endl;
	*my_file << "ASCII" << endl;
}

void output_vtk::write_grid(const grid* lgrid)
{
  int i;
  int nbNoe = lgrid->get_nbNoe();
  
  *my_file << "DATASET STRUCTURED_POINTS" << endl;
  
  *my_file << "DIMENSIONS ";
  for (i=0; i<3; i++)
    *my_file  << lgrid->get_nbNoe_xyz(i) << " ";
  *my_file << endl;
  
  *my_file << "ORIGIN 0 0 0" << endl;
  
  *my_file << "SPACING ";
  for (i=0; i<3; i++)
    *my_file  << lgrid->get_dxdydz(i) << " ";
  *my_file << endl;  
}

void output_vtk::write_scalarfield(int nbV, const char* name, const double * scalar_field)
{
  if ( scalar_field )
  {
    *my_file << "SCALARS "<< name << " FLOAT " << endl;
    *my_file << "LOOKUP_TABLE default "<< endl;
    for(int i=0;i<nbV; i++)
      *my_file << scalar_field[i] << endl;
    *my_file << endl;
  }
  else
  {
    ostringstream os;
    os << "Can not write scalar field " << name << " into the output file";
    parallel::print_message(os);
  }
}

void output_vtk::write_scalarfield(int nbV, const char* name, const int * scalar_field)
{
  if ( scalar_field )
  {
    *my_file << "SCALARS "<< name << " FLOAT " << endl;
    *my_file << "LOOKUP_TABLE default "<< endl;
    for(int i=0;i<nbV; i++)
      *my_file << scalar_field[i] << endl;
    *my_file << endl;
  }
    else
  {
    ostringstream os;
    os << "Can not write scalar field " << name << " into the output file";
    parallel::print_message(os);
  }
}

