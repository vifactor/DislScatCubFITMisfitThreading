/*
 * ProgramSettings.h
 *
 *  Created on: 11 june 2013
 *  Modified on: 17 june 2014
 *      Author: kopp
 */

#ifndef PROGRAMSETTINGS_H_
#define PROGRAMSETTINGS_H_

#include "StringTools.h"
#include "NonlinearFit.h"
#include <MillerIndexCub.h>
#include <libconfig.h++>
#include <boost/filesystem.hpp>

class ProgramSettings
{
public:
	class Exception: public std::exception
	{
	public:
		Exception(std::string m)
		{
			msg = "ProgramSettings::" + m;
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
	struct CalculatorSettings
	{
		/*cubic Miller indices*/
		MillerReciprocalCubIndices Q;
		/*X-ray wavelength*/
		double lambda;

		double qresolX, qresolZ;
		NonlinearFit::FitParameter scale;
		NonlinearFit::FitParameter background;

		long int sampling;
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
			double b_edge, b_screw;
			/*dislocation density*/
			double rho;
			/*correlation radius*/
			double rc;
		};

		MisfitDislocationType misfit;
		ThreadingDislocationType threading_edge;
		ThreadingDislocationType threading_screw;
		ThreadingDislocationType threading_mixed;
        friend std::ostream& operator<< (std::ostream &out,
                                        const SampleConfig &data);
        void set(const libconfig::Setting&, 
                NonlinearFit::CalculatorParameterMap&);
	};
	const CalculatorSettings& getCalculatorSettings() const
	{
		return m_calculatorSettings;
	}
	
	const SampleConfig& getSampleConfig() const
	{
		return m_sampleConfig;
	}
	const DataConfig& getDataConfig() const
	{
		return m_dataConfig;
	}
	const FitConfig& getFitConfig() const
	{
		return m_fitConfig;
	}
	const NonlinearFit::CalculatorParameterMap& getCPMap() const
	{
	    return m_cpMap;
	}
	const boost::filesystem::path& getConfigfile() const
	{
		return m_cfgfile;
	}
	const boost::filesystem::path& getResumefile() const
	{
		return m_resfile;
	}
	ProgramSettings();
	virtual ~ProgramSettings(){}

	void read(const boost::filesystem::path& cfgdir);
	void print() const;
protected:
	void readCalculatorSettings(const libconfig::Setting& root);

	void printCalculatorSettings() const;

	void printCoplanarSettings() const;


	CalculatorSettings m_calculatorSettings;

	boost::filesystem::path m_cfgfile;
	boost::filesystem::path m_resfile;
	
	NonlinearFit::CalculatorParameterMap m_cpMap;
	SampleConfig m_sampleConfig;
	DataConfig m_dataConfig;
	FitConfig m_fitConfig;
	
};

#endif /* PROGRAMSETTINGS_H_ */
