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
	ANACalculatorCoplanarTripleArgument(double qx = 0, double qz = 0) {m_qx = qx; m_qz = qz;}
	virtual ~ANACalculatorCoplanarTripleArgument() {}
	double m_qx, m_qz;
};

class FitANACalculatorCoplanarTriple:
		public NonlinearFit::FitCalculator
{
public:
	FitANACalculatorCoplanarTriple(ANACalculatorCoplanarTriple * calculator, ANASampleHex * sample);
	virtual ~FitANACalculatorCoplanarTriple();
	virtual void reinit(const NonlinearFit::CalculatorParameterMap& params);
	virtual double eval(const NonlinearFit::CalculatorArgument * arg);
protected:
	ANACalculatorCoplanarTriple * m_calculator;
	ANASampleHex * m_sample;
	double m_scale, m_background;
};

#endif /* FITANACALCULATORCOPLANARTRIPLE_H_ */
