//============================================================================
// Name        : DislScatANAMisfitHex.cpp
// Author      : Viktor Kopp
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include "Engine.h"

using namespace std;

void copy(std::string src, std::string dest);

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

void copy(std::string src, std::string dest)
{
	std::ofstream fout;
	std::ifstream fin;
	std::string filename, line;

	filename = dest;
	stripExtension(filename);
	filename += ".~cfg";

	fout.open(filename.c_str());
	fin.open(src.c_str());

	while (!fin.eof())
	{
		getline(fin, line);
		fout << line;
	}

	fin.close();
	fout.close();
}
