#ifndef _OUTPUT_H_
#define _OUTPUT_H_

#include <string>
#include <fstream>

#include "grid.h"

using namespace std;

class output
{
public:
	output();
	virtual ~output();

	int open_newfile();
	void close_file();

	void increment_counter() { counter++; }

	void set_filename(string & lname) { generic_name = lname; }
	void set_counter(int lcounteur) { counter = lcounteur; }
	void set_file(ofstream* lfile) { my_file = lfile; }

	const int& get_counter() const { return counter; }
	const string& get_fullname() const { return full_name; }
	const ofstream * get_file() const { return my_file; }

	virtual const string& get_extension() const { return extension; }
	virtual void write_grid(const grid* lgrid)=0;
	virtual void write_scalarfield(int nbV, const char* name, const double * scalar_field)=0;
	virtual void write_scalarfield(int nbV, const char* name, const int * scalar_field)=0;
	virtual void write_header_point_data(int nbV) {;}

protected:
	string generic_name;
	string full_name;
	static const string extension;
	ofstream* my_file;
	int counter;

	virtual void write_header() { ; }	
};

#endif
