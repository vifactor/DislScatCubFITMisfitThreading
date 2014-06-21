/*
 * FitANACalculatorCoplanarTriple.cpp
 *
 *  Created on: 25 july 2013
 *      Author: kopp
 */

#include "FitANACalculatorCoplanarTriple.h"
using namespace boost;

FitANACalculatorCoplanarTriple::FitANACalculatorCoplanarTriple(
                    const std::vector<ANACalculatorCoplanarTriple *>& calculators,
                    ANASampleCub * sample)
{
	m_sample = sample;
	
	m_calculators = calculators;
	m_backgrounds.resize(m_calculators.size());
    m_scales.resize(m_calculators.size());
    initParameterNames();
}

void FitANACalculatorCoplanarTriple::initParameterNames()
{
    m_scale_names.resize(m_calculators.size());
    m_background_names.resize(m_calculators.size());
    for(size_t id = 0; id < m_calculators.size(); ++id)
    {
        m_scale_names[id] = "Data.[" + lexical_cast<std::string>(id) + "].I0";
        m_background_names[id] = "Data.[" + lexical_cast<std::string>(id) + "].Ibg";
    }
}

void FitANACalculatorCoplanarTriple::reinit(const NonlinearFit::CalculatorParameterMap& params)
{
	double rho_mf;
	double rho_th_edge, rho_th_screw, rho_th_mixed;
	double rc_th_edge, rc_th_screw, rc_th_mixed;

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

	/*std::cout << "rho_th_edge:\t" << rho_th_edge << std::endl;
	std::cout << "rho_th_screw:\t" << rho_th_screw << std::endl;
	std::cout << "rho_th_mixed:\t" << rho_th_mixed << std::endl;

	std::cout << "rc_th_edge:\t" << rc_th_edge << std::endl;
	std::cout << "rc_th_screw:\t" << rc_th_screw << std::endl;
	std::cout << "rc_th_mixed:\t" << rc_th_mixed << std::endl;*/

	m_sample->resetThreadingLayer(0, rho_th_edge, rc_th_edge);
	m_sample->resetThreadingLayer(1, rho_th_screw, rc_th_screw);
	m_sample->resetThreadingLayer(2, rho_th_mixed, rc_th_mixed);
    
    for(size_t id = 0; id < m_calculators.size(); ++id)
    {
        /*reinitialization scale and background coefficients*/
        m_scales[id] = params.find(m_scale_names[id])->second;
        m_backgrounds[id] = params.find(m_background_names[id])->second;
        
        std::cout << m_scale_names[id] << ":\t" << m_scales[id] << std::endl;
        
        /*update each calculator*/
        m_calculators[id]->setSample(m_sample);
	}
}

double FitANACalculatorCoplanarTriple::eval(const NonlinearFit::CalculatorArgument * arg)
{
	static double qx, qz, result;
	static int id;

	qx = static_cast<const ANACalculatorCoplanarTripleArgument* >(arg)->m_qx;
	qz = static_cast<const ANACalculatorCoplanarTripleArgument* >(arg)->m_qz;
	id = static_cast<const ANACalculatorCoplanarTripleArgument* >(arg)->m_id;

	result = m_scales[id] * m_calculators[id]->I(qx, qz) + m_backgrounds[id];

	return result;
}

