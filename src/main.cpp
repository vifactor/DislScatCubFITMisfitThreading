//============================================================================
// Name        : DislScatANAMisfitHex.cpp
// Author      : Viktor Kopp
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "Engine.h"

using namespace std;

int main()
{
	Engine engine;
	std::string cfgfile;
	try
	{
		cfgfile = "default.cfg";
		engine.exec(cfgfile);
	}catch(const ProgramSettings::Exception& ex)
	{
		std::cout << ex.what() << std::endl;
		return -1;
	}

	std::cout << "Done." << std::endl;

	return 0;
}
