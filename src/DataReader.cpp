/*
 * DataReader.cpp
 *
 *  Created on: 5 mar. 2012
 *  Modified on: 21 june 2014
 *      Author: Viktor Kopp
 */

#include "DataReader.h"

DataReader::DataReader(const char * sep, char com)
{
	m_cs = new boost::char_separator<char>(sep);
	m_com = com;
}

DataReader::~DataReader()
{
    if(m_cs) delete m_cs;
}

void DataReader::init()
{
    m_data.clear();
    m_columnNames.clear();
}

void DataReader::parse(const std::string& filename)
{
    std::ifstream fin;
	std::string line;
	std::vector<std::string> values;
    
    /*clean containers*/
    init();
    
    /*throw exception if opening/reading/closing file failed*/
    fin.open(filename.c_str());
    if(!fin.is_open()) throw std::ifstream::failure(filename + " does not exists.");
    
    assignHeader(fin);
	while (getline(fin, line))
	{
		//split line into columns
		Tokenizer tok(line, *m_cs);
		values.assign(tok.begin(),tok.end());
		//check if nb of columns corresponds to number of headers
		if(values.size() == m_columnNames.size())
		{
		    //distribute everything by columns
            for(size_t i = 0; i < values.size(); ++i)
            {
                m_data[m_columnNames[i]].push_back(
                        boost::lexical_cast<double>(values[i]));
            }
        }
	}
	fin.close();
}

bool DataReader::isComment(const std::string& str) const
{
	//empty string
	if(str.empty()) return true;
	//comment
	if(str[0] == m_com) return true;
	return false;
}

void DataReader::assignHeader(std::ifstream & fin)
{
        /*returns a commented string preceeding uncommented ones*/
    std::string header, line;
    int pos;
    
    /*memorize initial position*/
    pos = fin.tellg();
    while(getline(fin,line))
    {
        /*remove spaces from front and back of the string*/
        boost::algorithm::trim(line);
        if (isComment(line))
        {
            header = line.substr(1);
            /*memorize initial position*/
            pos = fin.tellg();
        }
        else
        {
            fin.seekg(pos);
            break;
        }
    }
    
    Tokenizer tok(header, *m_cs);
    //split header line into column names
    m_columnNames.assign(tok.begin(), tok.end());
}

const DataReader::ColumnType& DataReader::columnGet(const ColumnName& name) const
{
	std::map<ColumnName, ColumnType>::const_iterator it_map;

	it_map = m_data.find(name);
	if(it_map != m_data.end())
	{
        return it_map->second;
	}
	else
	    //empty column
	    return m_emptyColumn;
}

bool DataReader::columnExists(const ColumnName& name) const
{
	std::map<ColumnName, ColumnType>::const_iterator it_map;

	it_map = m_data.find(name);
	if(it_map != m_data.end())
		return true;
	else
		return false;
}
