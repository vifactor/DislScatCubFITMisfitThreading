/*
 * FitANACalculatorCoplanarTriple.cpp
 *
 *  Created on: 25 лип. 2013
 *      Author: kopp
 */

#include "FitANACalculatorCoplanarTriple.h"

FitANACalculatorCoplanarTriple::FitANACalculatorCoplanarTriple(ANACalculatorCoplanarTriple * calculator, ANASampleHex * sample)
{
	m_calculator = calculator;
	m_sample = sample;
	m_background = 0.0;
	m_scale = 1.0;
}

FitANACalculatorCoplanarTriple::~FitANACalculatorCoplanarTriple()
{
}

void FitANACalculatorCoplanarTriple::reinit(const NonlinearFit::CalculatorParameterMap& params)
{
	double rho_mf;
	double rho_th_edge, rho_th_screw, rho_th_mixed;
	double rc_th_edge, rc_th_screw, rc_th_mixed;

	/*reinitialization of correlation radii of threading dislocations*/
	m_scale = params.find("Calculator.scale")->second;
	m_background = params.find("Calculator.background")->second;

	/*
	 * reinitialization of densities of misfit dislocations
	 * initially threading dislocation density is given in [cm-1]
	 * coefficient 1e-7 transforms it to [nm-1]
	*/
	rho_mf = params.find("Sample.dislocations.misfit.rho")->second  * 1e-7;

	m_sample->resetMisfitInterface(0, rho_mf);

	/*
	 * reinitialization of densities of threading dislocations
	 * initially threading dislocation density is given in [cm-2]
	 * coefficient 1e-14 transforms it to [nm-2]
	*/
	rho_th_edge = params.find("Sample.dislocations.threading.edge.rho")->second  * 1e-14;
	rho_th_screw = params.find("Sample.dislocations.threading.screw.rho")->second * 1e-14;
	rho_th_mixed = params.find("Sample.dislocations.threading.mixed.rho")->second * 1e-14;

	/*reinitialization of correlation radii of threading dislocations*/
	rc_th_edge = params.find("Sample.dislocations.threading.edge.rc")->second;
	rc_th_screw = params.find("Sample.dislocations.threading.screw.rc")->second;
	rc_th_mixed = params.find("Sample.dislocations.threading.mixed.rc")->second;

	std::cout << "scale:\t" << m_scale << std::endl;
	std::cout << "rho_mf:\t" << rho_mf << std::endl;

	/*std::cout << "rho_th_edge:\t" << rho_th_edge << std::endl;
	std::cout << "rho_th_screw:\t" << rho_th_screw << std::endl;
	std::cout << "rho_th_mixed:\t" << rho_th_mixed << std::endl;

	std::cout << "rc_th_edge:\t" << rc_th_edge << std::endl;
	std::cout << "rc_th_screw:\t" << rc_th_screw << std::endl;
	std::cout << "rc_th_mixed:\t" << rc_th_mixed << std::endl;*/

	m_sample->resetThreadingLayer(0, rho_th_edge, rc_th_edge);
	m_sample->resetThreadingLayer(1, rho_th_screw, rc_th_screw);
	m_sample->resetThreadingLayer(2, rho_th_mixed, rc_th_mixed);

	m_calculator->setSample(m_sample);
}

double FitANACalculatorCoplanarTriple::eval(const NonlinearFit::CalculatorArgument * arg)
{
	static double qx, qz, result;

	qx = static_cast<const ANACalculatorCoplanarTripleArgument* >(arg)->m_qx;
	qz = static_cast<const ANACalculatorCoplanarTripleArgument* >(arg)->m_qz;

	result = m_scale * m_calculator->I(qx, qz) + m_background;

	return result;
}

