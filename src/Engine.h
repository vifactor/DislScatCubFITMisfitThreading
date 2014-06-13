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
		Exception(std::string m) :
			msg("Engine::" + m){}
		~Exception() throw () {}
		const char* what() const throw () {return msg.c_str();}
	private:
		std::string msg;
	};
	Engine();
	void exec(const std::string & cfgname);
	virtual ~Engine();
private:
	void init();
	void setupComponents();
	void setupCalculator();
	void setupCParameters();
	void setupFitter();
	void addFParameter(const NonlinearFit::FitParameter & param);
	void readData();
	void fillData();

	void doWork();
	void calcI(std::vector<double>& pts);
	void fitI();

	std::string getFilename() const;
	void saveSettings() const;
	void saveResume() const;
	void saveResult() const;
	void saveCalcIntensity() const;
	void saveFitIntensity() const;
	void printFitterInfo(NonlinearFit::NonlinearFitter * fitter);

	ANASampleHex  * m_sample;
	ANACalculatorCoplanarTriple *  m_calculator;
	FitANACalculatorCoplanarTriple *  m_fit_calculator;
	NonlinearFit::NonlinearFitter * m_fitter;
	ProgramSettings * m_programSettings;

	/*vector of fit parameters with their initial values*/
	NonlinearFit::FitParameterList m_fParametersInit;
	/*vector of fit parameters with their final values*/
	NonlinearFit::FitParameterList m_fParametersFinal;
	/*map of calculator parameters (potential fit parameters)*/
	NonlinearFit::CalculatorParameterMap m_cParameters;

	std::vector<double> m_qx_vals, m_qz_vals;
	std::vector<double> m_exp_intens_vals, m_ini_intens_vals, m_fin_intens_vals;
	NonlinearFit::DataPointList m_DataPoints;
	NonlinearFit::ResidualList m_Residuals;
};

#endif /* Engine_H_ */
