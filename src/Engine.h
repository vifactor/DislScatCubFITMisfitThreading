/*
 * Engine.h
 *
 *  Created on: 30 may 2012
 *      Author: kopp
 */

#ifndef Engine_H_
#define Engine_H_

#include "FitANACalculatorCoplanarTriple.h"
#include "PortFitter.h"
#include "ProgramSettings.h"
#include "DataReader.h"

#include <fstream>

class Engine
{
public:
	class Exception: public std::exception
	{
	public:
		Exception(const std::string& m) :
			msg("Engine::" + m){}
		~Exception() throw () {}
		const char* what() const throw () {return msg.c_str();}
	private:
		std::string msg;
	};
	Engine();
	void exec(const boost::filesystem::path & cfgfolder);
	virtual ~Engine();
private:
	void init();
	void setupComponents();
	void setupCalculator();
	void setupFitter();
	void readData();

	void doWork();
	void calcI(std::vector<double>& pts);
	void fitI();

	boost::filesystem::path getOutFilename(
	                const boost::filesystem::path & filename) const;
	void updateConfigFile(const boost::filesystem::path & filename) const;
	void saveResume() const;
	void saveResult() const;
	void saveFitIntensity() const;
	void printFitterInfo(NonlinearFit::NonlinearFitter * fitter);

	ANASampleCub  * m_sample;
	std::vector<ANACalculatorCoplanarTriple *>  m_calculators;
	FitANACalculatorCoplanarTriple *  m_fit_calculator;
	NonlinearFit::NonlinearFitter * m_fitter;
	ProgramSettings * m_programSettings;

	/*vector of fit parameters with their final values*/
	NonlinearFit::FitParameterList m_fParametersFinal;

	std::vector<double> m_qx_vals, m_qz_vals;
	std::vector<double> m_exp_intens_vals, m_ini_intens_vals, m_fin_intens_vals;
	NonlinearFit::DataPointList m_DataPoints;
	NonlinearFit::ResidualList m_Residuals;
	
	boost::filesystem::path m_WorkDir;
};

#endif /* Engine_H_ */
