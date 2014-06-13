/*
 * StringTools.h
 *
 *  Created on: 12 feb. 2009
 *  Modified on: 13 aug. 2012
 *	Modified on: 10 dec. 2012
 *      Author: Viktor Kopp
 */

#ifndef STRINGTOOLS_H_
#define STRINGTOOLS_H_

#include<string>
#include<vector>
#include<sstream>
#include <cctype>
#include <iostream>
#include <algorithm>

//splits string into two parts:
//returns part after the last point (extension)
//argument transformed into part before last point(base name)
extern std::string stripExtension(std::string& filename);
extern std::size_t stripFirstBlanc(std::string& str);
extern std::size_t stripAllBlanc(std::string& str);

//Splits string "s" with a list of delimiters in "delims"
extern std::vector<std::string> split(std::string s, std::string delims);

/*very simple parser utilities for command line arguments*/
extern char* getCmdOption(char ** begin, char ** end, const std::string & option);
extern bool cmdOptionExists(char** begin, char** end, const std::string& option);

//transforms any object which has a "<<" operator into "std::string"
template<typename T>
std::string toString(T t) 
{
    std::stringstream s;
    s << t;
    return s.str();
}

#endif /* STRINGTOOLS_H_ */
