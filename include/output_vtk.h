
#ifndef _OUTPUT_VTK_H_
#define _OUTPUT_VTK_H_

#include "output.h"

class output_vtk : public output
{
public:

	output_vtk() : output() { ; }

	virtual const string& get_extension() const;
	virtual void write_grid(const grid* lgrid);
	virtual void write_scalarfield(int nbV, const char* name, const double * scalar_field);
	virtual void write_scalarfield(int nbV, const char* name, const int * scalar_field);
	virtual void write_header_point_data(int nbV);

protected:

	virtual void write_header();
	static const string extension;
};


#endif
