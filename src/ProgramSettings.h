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
	struct SampleSettings
	{
		double nu;
		double thickness;
		double width;
		/*hexagonal lattice parameters*/
		double a0;

		struct MisfitDislocationType
		{
			/*burgers components*/
			MillerDirectCubIndices b, l;
			NonlinearFit::FitParameter rho;
		};
		struct ThreadingDislocationType
		{
			double b_edge, b_screw;
			NonlinearFit::FitParameter rho;
			NonlinearFit::FitParameter rc;
		};

		MisfitDislocationType misfit;
		ThreadingDislocationType threading_edge;
		ThreadingDislocationType threading_screw;
		ThreadingDislocationType threading_mixed;
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
	struct EngineSettings
	{
		std::string outfile;
		std::string datafile;
		std::string resumefile;
		int nb_iter;
	};
	struct DataConfig
	{
        MillerReciprocalCubIndices Q;
        double resolX, resolZ;
        boost::filesystem::path file;
        double Ibg, I0;
        
        friend std::ostream& operator<< (std::ostream &out, const DataConfig &data);
        void set(const libconfig::Setting&);
	};
	const SampleSettings& getSampleSettings() const
	{
		return m_sampleSettings;
	}
	const CalculatorSettings& getCalculatorSettings() const
	{
		return m_calculatorSettings;
	}
	const EngineSettings& getEngineSettings() const
	{
		return m_engineSettings;
	}
	const boost::filesystem::path& getConfigfile() const
	{
		return m_cfgfile;
	}
	ProgramSettings();
	virtual ~ProgramSettings();

	void read(const boost::filesystem::path& cfgdir);
	void print() const;
protected:
	void readCalculatorSettings(const libconfig::Setting& root);
	void readSampleSettings(const libconfig::Setting& root);
	void readEngineSettings(const libconfig::Setting& root);

	void readMisfitDislocations(const libconfig::Setting& stg);
	void readThreadingDislocations(const libconfig::Setting& stg);

	void printSampleSettings() const;
	void printCalculatorSettings() const;
	void printEngineSettings() const;

	void printCoplanarSettings() const;

	void printMisfitDislocations() const;
	void printThreadingDislocations() const;

	SampleSettings m_sampleSettings;
	CalculatorSettings m_calculatorSettings;
	EngineSettings m_engineSettings;
	boost::filesystem::path m_cfgfile;
};

#endif /* PROGRAMSETTINGS_H_ */
