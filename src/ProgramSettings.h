/*
 * ProgramSettings.h
 *
 *  Created on: 11 june 2013
 *  Modified on: 17 june 2014
 *      Author: kopp
 */

#ifndef PROGRAMSETTINGS_H_
#define PROGRAMSETTINGS_H_

#include "NonlinearFit.h"
#include <MillerIndexCub.h>
#include <libconfig.h++>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

class ProgramSettings
{
public:
	class Exception: public std::exception
	{
	public:
		Exception(const std::string& m)
		{
			msg = "ProgramSettings::" + m;
		}
		Exception(const char * m)
		{
			msg = "ProgramSettings::";
			msg += m;
		}
		~Exception() throw ()
		{
		}
		const char* what() const throw ()
		{
			return msg.c_str();
		}
	private:
		std::string msg;
	};

	struct DataConfig
	{
        MillerReciprocalCubIndices Q;
        double resolX, resolZ;
        boost::filesystem::path file;
        double Ibg, I0;
        
        friend std::ostream& operator<< (std::ostream &out, const DataConfig &data);
        void set(const libconfig::Setting&, 
                NonlinearFit::CalculatorParameterMap&);
	};
	struct FitConfig
	{
	    int nbIterations;
	    NonlinearFit::FitParameterList fitParameters;
	    
        friend std::ostream& operator<< (std::ostream &out,
                                        const FitConfig &data);
        void set(const libconfig::Setting&, 
                const NonlinearFit::CalculatorParameterMap&);
	};
	struct SampleConfig
	{
	    /*Poisson ratio*/
	    double nu;
		/*sample dimensions*/
		double thickness;
		double width;
		/*cubic lattice parameters*/
		double a0;

		struct MisfitDislocationType
		{
			/*burgers components*/
			MillerDirectCubIndices b, l;
			double rho;
		};
		struct ThreadingDislocationType
		{
			MillerDirectCubIndices b;
			/*dislocation density*/
			double rho;
			/*correlation radius*/
			double rc;
		};

        std::vector<MisfitDislocationType> misfit;
		std::vector<ThreadingDislocationType> threading;
        friend std::ostream& operator<< (std::ostream &out,
                                        const SampleConfig &data);
        void set(const libconfig::Setting&, 
                NonlinearFit::CalculatorParameterMap&);
	};
	
	const SampleConfig& getSampleConfig() const
	{
		return m_sampleConfig;
	}
	const std::vector<DataConfig>& getDataConfig() const
	{
		return m_dataConfig;
	}
	const DataConfig& getDataConfig(size_t i) const
	{
		return m_dataConfig[i];
	}
	const FitConfig& getFitConfig() const
	{
		return m_fitConfig;
	}
	const NonlinearFit::CalculatorParameterMap& getCPMap() const
	{
	    return m_cpMap;
	}
	NonlinearFit::CalculatorParameterMap& getCPMap()
	{
	    return m_cpMap;
	}
	const boost::filesystem::path& getDataConfigFile() const
	{
		return m_datfile;
	}
		const boost::filesystem::path& getFitConfigFile() const
	{
		return m_fitfile;
	}
	const boost::filesystem::path& getSampleConfigFile() const
	{
		return m_samfile;
	}
	const boost::filesystem::path& getResumefile() const
	{
		return m_resfile;
	}
	ProgramSettings();
	virtual ~ProgramSettings(){}

	void read(const boost::filesystem::path& cfgdir);
protected:

	boost::filesystem::path m_resfile;
	boost::filesystem::path m_datfile;
	boost::filesystem::path m_samfile;
	boost::filesystem::path m_fitfile;
	
	NonlinearFit::CalculatorParameterMap m_cpMap;
	SampleConfig m_sampleConfig;
	std::vector<DataConfig> m_dataConfig;
	FitConfig m_fitConfig;
	
};

#endif /* PROGRAMSETTINGS_H_ */
