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
	ANACalculatorCoplanarTripleArgument(double qx, double qz,
	                                    ANACalculatorCoplanarTriple * calculator)
	{
	    m_qx = qx; m_qz = qz;
	    m_calculator = calculator;
	}
	virtual ~ANACalculatorCoplanarTripleArgument() {}
	
	double m_qx, m_qz;
	ANACalculatorCoplanarTriple * m_calculator;
};

class FitANACalculatorCoplanarTriple:
		public NonlinearFit::FitCalculator
{
public:
	FitANACalculatorCoplanarTriple(ANACalculatorCoplanarTriple * calculator, ANASampleCub * sample);
	virtual ~FitANACalculatorCoplanarTriple() {}
	virtual void reinit(const NonlinearFit::CalculatorParameterMap& params);
	virtual double eval(const NonlinearFit::CalculatorArgument * arg);
protected:
	ANACalculatorCoplanarTriple * m_calculator;
	ANASampleCub * m_sample;
	double m_scale, m_background;
};

#endif /* FITANACALCULATORCOPLANARTRIPLE_H_ */
