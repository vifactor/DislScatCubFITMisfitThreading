/*
 * StringTools.cpp
 *
 *  Created on: 12 feb. 2009
 *  Modified on: 13 aug. 2012
 *	Modified on: 10 dec. 2012
 *      Author: Viktor Kopp
 */

#include "StringTools.h"

std::string stripExtension(std::string& filename)
{
	std::string ext="";//no extension
	// find last '.' in the file
	std::string::size_type pos(filename.find_last_of('.'));
	if (pos != filename.npos)//has an extension
    {
		ext.assign(filename, pos+1, filename.npos);
		// put period into extension
		filename.assign(filename, 0, pos);
    }
	return ext;
}

std::size_t stripFirstBlanc(std::string& str)
{
	std::size_t nzeros=0;
	while(isspace(str[nzeros++]));
	str.erase(0, nzeros-1);
	return nzeros;
}

std::size_t stripAllBlanc(std::string& str)
{
	std::size_t nzeros=0;
	std::string temp("");
	for(std::size_t i=0; i<str.length(); i++)
	{
		if(isspace(str[i]))
	    {
			nzeros++;
	    }
		else
		{
			temp+=str[i];
		}
	}
	str=temp;
	return nzeros;
}

/* ++ 13 august 2012 */

//Splits string "s" with a list of delimiters in "delims"
std::vector<std::string> split(std::string str, std::string delims)
{
	std::vector<std::string> result;
	size_t lastPos = 0;
	size_t pos = str.find_first_of(delims);

	while (pos != std::string::npos)
	{
		if (pos != lastPos)
			result.push_back(str.substr(lastPos, pos - lastPos));
		lastPos = pos + 1;
		pos = str.find_first_of(delims, lastPos);
	}
	if (lastPos < str.length())
		result.push_back(str.substr(lastPos, pos - lastPos));

	return result;
}

/*++ 10 december 2012*/
//command line utilities
char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}
