
#include <sstream>
#include <iostream>

#include "output.h"

const string output::extension = "";


output::output()
{
	counter = 0;
	my_file = 0;
	generic_name = "";
	full_name = "";
}

output::~output()
{
	delete my_file;
}

int output::open_newfile()
{
	int ok = 1;

	if ( my_file )
		close_file();

	my_file = new ofstream();
	if ( my_file )
	{
		ostringstream lfull_name;
		lfull_name << generic_name << "_" << counter << "." << get_extension();
		full_name = lfull_name.str();
		my_file->open(full_name.c_str(), ios::out);
		if ( my_file->is_open() )
			write_header();
		else
			ok = 0;
	}
	else
		ok = 0;

	return ok;
}

void output::close_file()
{
	if ( my_file )
	{
		my_file->close();
		delete my_file;
		my_file = 0;
	}
}
