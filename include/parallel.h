
#ifndef _PARALLEL_H_
#define _PARALLEL_H_

#include <mpi.h>
#include <stdlib.h>
#include <sstream>


namespace parallel
{
	const int master = 0;
	const MPI_Comm comm = MPI_COMM_WORLD;

	inline int my_rank()
	{
		int rank;
		MPI_Comm_rank(comm,&rank);
		return rank;
	}

	inline int nb_proc()
	{
		int size;
		MPI_Comm_size(comm,&size);
		return size;
	}
	
	inline bool IamMaster()
	{
	 return (my_rank() == master); 
	}

	inline void print_message(const char* os)
	{
		if ( my_rank() == master )
			std::cout << os << std::endl;
	}

	inline void print_message(std::ostringstream & os)
	{
		if ( my_rank() == master )
			std::cout << os.str().c_str() << std::endl;
	}

	inline void print_problem(const char* os)
	{
		if ( my_rank() == master )
			std::cout << "Problem! " << os << std::endl;
		exit(EXIT_FAILURE);
	}

	inline void print_problem(std::ostringstream & os)
	{
		if ( my_rank() == master )
			std::cout << "Problem! " << os.str().c_str() << std::endl;
		exit(EXIT_FAILURE);
	}

}

#endif
