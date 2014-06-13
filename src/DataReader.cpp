/*
 * DataReader.cpp
 *
 *  Created on: 5 mar. 2012
 *      Author: kopp
 */

#include "DataReader.h"

DataReader::DataReader()
{
	nbRows = 0;
	status = true;
}

DataReader::~DataReader()
{
}

void DataReader::readFile(std::string filename)
{
	this->filename=filename;

	std::ifstream fin(filename.c_str());
	if (!fin)
	{
		status = false;
		return;
	}

	std::string line;
	std::istringstream is;
	size_t ln;
	std::vector<std::string> values;
	double value;

	//here we read a header line
	getline(fin, line);
	readHeader(line);

	ln = 1;
	while (!fin.eof())
	{
		getline(fin, line);
		++ln;
		//skip the comment
		if (isComment(line))
			continue;

		//split line into columns
		values = split(line, " \t\n");
		//check if nb of columns corresponds to number of headers
		if(values.size() != columnNames.size())
		{
			continue;
		}

		for(size_t icol = 0; icol < columnNames.size(); icol++)
		{
			is.str(values[icol]);
			is >> value;
			is.clear();

			data[columnNames[icol]].push_back(value);
			dataf[columnNames[icol]].push_back(value);

		}
		++nbRows;

	}
	fin.close();
}

bool DataReader::isComment(const std::string& str)
{
	//empty string
	if(str.empty()) return true;
	//comment
	if(str[0]=='#') return true;
	return false;
}

void DataReader::readHeader(const std::string& str)
{
	//no header line of type: % colName1 colName2...
	if(str[0] != '#')
	{
		status = false;
		return;
	}

	//split header line into column names with delimiters "% \t\n"
	columnNames = split(str, "# \t\n\r");

	if(columnNames.empty())
	{
		status = false;
		return;
	}
}

void DataReader::getColumn(ColumnType& col, const ColumnName& name) const
{
	std::map<ColumnName, ColumnType>::const_iterator it_map;

	col.clear();

	it_map = dataf.find(name);
	if(it_map != dataf.end())
	{
		for(size_t i = 0; i < it_map->second.size(); ++i)
			col.push_back(it_map->second.at(i));
	}
	else
	{
		status = false;
	}
}

bool DataReader::columnExist(const ColumnName& name) const
{
	std::map<ColumnName, ColumnType>::const_iterator it_map;

	it_map = dataf.find(name);
	if(it_map != dataf.end())
	{
		return true;
	}
	else
	{
		return false;
	}
}
