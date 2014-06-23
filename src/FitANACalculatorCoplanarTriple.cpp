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
	m_scales.resize(m_calculators.size());
    m_backgrounds.resize(m_calculators.size());
    
    initParameterNames();
}

void FitANACalculatorCoplanarTriple::initParameterNames()
{
    /*Data fit parameter names*/
    m_scale_names.resize(m_calculators.size());
    m_background_names.resize(m_calculators.size());
    for(size_t id = 0; id < m_calculators.size(); ++id)
    {
        m_scale_names.at(id) = "Data.[" + lexical_cast<std::string>(id) + "].I0";
        m_background_names.at(id) = "Data.[" + lexical_cast<std::string>(id) + "].Ibg";
    }
    /*misfit interfaces names*/
    m_mf_density_names.resize(m_sample->getNbMisfitInterfaces());
    for(size_t id = 0; id < m_sample->getNbMisfitInterfaces(); ++id)
    {
        m_mf_density_names.at(id) = "Sample.dislocations.misfit.[" 
                            + lexical_cast<std::string>(id) 
                            + "].rho";
    }
    /*theading layers names*/
    m_th_density_names.resize(m_sample->getNbThreadingLayers());
    m_th_rc_names.resize(m_sample->getNbThreadingLayers());
    for(size_t id = 0; id < m_sample->getNbThreadingLayers(); ++id)
    {
        m_th_density_names.at(id) = "Sample.dislocations.threading.[" 
                            + lexical_cast<std::string>(id) 
                            + "].rho";
        m_th_rc_names.at(id) = "Sample.dislocations.threading.[" 
                            + lexical_cast<std::string>(id) 
                            + "].rc";
    }
}

void FitANACalculatorCoplanarTriple::reinit(const NonlinearFit::CalculatorParameterMap& params)
{
    static double rho_mf;
	static double rho_th, rc_th;

	/*
	 * reinitialization of densities of misfit dislocations.
	 * initially misfit dislocation density is given in [cm-1]
	 * coefficient 1e-7 transforms it to [nm-1]
	*/
	for(size_t id = 0; id < m_sample->getNbMisfitInterfaces(); ++id)
	{
	    rho_mf = params.find(m_mf_density_names[id])->second  * 1e-7;
	    m_sample->resetMisfitInterface(id, rho_mf);
	}

	/*
	 * reinitialization of densities and correlation radii of threading dislocations.
	 * initially threading dislocation density is given in [cm-2]
	 * coefficient 1e-14 transforms it to [nm-2]
	*/
	for(size_t id = 0; id < m_sample->getNbThreadingLayers(); ++id)
	{
	    rho_th = params.find(m_th_density_names[id])->second * 1e-14;
	    rc_th = params.find(m_th_rc_names[id])->second;
	    m_sample->resetThreadingLayer(id, rho_th, rc_th);
	    
	    /*//DEBUG*/
	    std::cout << m_th_density_names[id] << ":\t" << rho_th << std::endl;
	    std::cout << m_th_rc_names[id]  << ":\t" << rc_th << std::endl;
	}
    
    for(size_t id = 0; id < m_calculators.size(); ++id)
    {
        /*reinitialization scale and background coefficients*/
        m_scales[id] = params.find(m_scale_names[id])->second;
        m_backgrounds[id] = params.find(m_background_names[id])->second;
        
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

