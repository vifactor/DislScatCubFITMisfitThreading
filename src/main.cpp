//============================================================================
// Name        : DislScatANAMisfitHex.cpp
// Author      : Viktor Kopp
// Version     :
// Copyright   : Your copyright notice
// Description : Fit of X-ray data to model including 
//              misfit and threading dislocation
//============================================================================

#include "Engine.h"

int main(int argc, char ** argv)
{
	if(argc > 1)
	{
	    Engine engine;
	    try
	    {
	        /* argv[1] contains name of the work folder*/
		    engine.exec(argv[1]);
	    }catch(const ProgramSettings::Exception& ex)
	    {
		    std::cout << ex.what() << std::endl;
		    return EXIT_FAILURE;
	    }
	}
	else
	{
	    std::cout << "No work folder specified." << std::endl;
	}

	std::cout << "Done." << std::endl;

	return EXIT_SUCCESS;
}
