/*
 * FitANACalculatorCoplanarTriple.h
 *
 *  Created on: 25 july 2013
 *      Author: kopp
 */

#ifndef FITANACALCULATORCOPLANARTRIPLE_H_
#define FITANACALCULATORCOPLANARTRIPLE_H_

#include "ANACalculatorCoplanarTriple.h"
#include "NonlinearFit.h"

class ANACalculatorCoplanarTripleArgument : public NonlinearFit::CalculatorArgument
{
public:
	ANACalculatorCoplanarTripleArgument(double qx, double qz, int id)
	{
	    m_qx = qx; m_qz = qz;
	    m_id = id;
	}
	virtual ~ANACalculatorCoplanarTripleArgument() {}
	
	double m_qx, m_qz;
	int m_id;
};

class FitANACalculatorCoplanarTriple:
		public NonlinearFit::FitCalculator
{
public:
	FitANACalculatorCoplanarTriple(
	        const std::vector<ANACalculatorCoplanarTriple *>& calculators,
	        ANASampleCub * sample);
	virtual ~FitANACalculatorCoplanarTriple() {}
	virtual void reinit(const NonlinearFit::CalculatorParameterMap& params);
	virtual double eval(const NonlinearFit::CalculatorArgument * arg);
protected:
    ANASampleCub * m_sample;
	std::vector<ANACalculatorCoplanarTriple *> m_calculators;
	std::vector<double> m_scales, m_backgrounds;
};

#endif /* FITANACALCULATORCOPLANARTRIPLE_H_ */
