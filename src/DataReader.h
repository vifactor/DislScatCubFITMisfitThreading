/*
 * DataReader.h
 *
 *  Created on: 5 mar. 2012
 *      Author: kopp
 */

#ifndef DATAREADER_H_
#define DATAREADER_H_

#include <fstream>
#include <map>
#include <algorithm>

#include <boost/algorithm/string/trim.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

class DataReader
{
public:
	typedef std::string ColumnName;
	typedef std::vector<double> ColumnType;
	typedef std::map<ColumnName, double> RowType;
	typedef boost::tokenizer<boost::char_separator<char> > Tokenizer;
public:
	DataReader(const char * sep = " ", char com = '#');
	virtual ~DataReader();
    void parse(const std::string& filename);
	bool columnExists(const ColumnName& name) const;
	const DataReader::ColumnType& columnGet(const ColumnName& name) const;
private:
    std::map<ColumnName, ColumnType> m_data;
    std::vector<ColumnName> m_columnNames;
    ColumnType m_emptyColumn;
    boost::char_separator<char> * m_cs;
    char m_com;

    bool isComment(const std::string& str) const;
    void assignHeader(std::ifstream & fin);
    void init();
};

#endif /* DATAREADER_H_ */
