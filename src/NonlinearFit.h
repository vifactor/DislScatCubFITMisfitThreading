/*
 * NonlinearFitter.h
 *
 *  Created on: 12 july 2012
 *      Author: Viktor Kopp
 *
 *      Provides a convenient interface to various nonlinear fit routines
 *      with simple boundary conditions
 *
 *      Version 0.3 from 120723
 */

#ifndef NONLINEARFIT_H_
#define NONLINEARFIT_H_

#include <map>
#include <vector>
#include <string>
#include <memory>

#include <iostream>

#include <gsl/gsl_math.h>

namespace NonlinearFit
{

/*name of a parameter is a string*/
typedef std::string NameType;
/*parameters are of type double*/
typedef double ParameterValueType;
typedef double ResidualValueType;
typedef std::vector<ResidualValueType> ResidualList;
typedef std::map<NameType, ParameterValueType> CalculatorParameterMap;
class CalculatorArgument
{
public:
	CalculatorArgument() {}
	virtual ~CalculatorArgument() {}
private:
};
struct DataPoint
{
	DataPoint(CalculatorArgument * arg, ResidualValueType resid) :
		m_Argument(arg), m_Residual(resid) {}
	CalculatorArgument * m_Argument;
	ResidualValueType m_Residual;
};
typedef std::vector<DataPoint> DataPointList;
struct FitParameter
{
	FitParameter(NameType name="", ParameterValueType value = 0.0,
			ParameterValueType lbvalue = 0.0, ParameterValueType ubvalue = 0.0,
			ParameterValueType scvalue = 1.0) :
			m_Name(name), m_Value(value), m_Lbvalue(lbvalue), m_Ubvalue(
					ubvalue), m_Scvalue(scvalue){}
	NameType m_Name;
	ParameterValueType m_Value;
	ParameterValueType m_Lbvalue;
	ParameterValueType m_Ubvalue;
	ParameterValueType m_Scvalue;
};

typedef std::vector<FitParameter> FitParameterList;

class FitCalculator
{
public:
	FitCalculator(){}
	virtual ~FitCalculator(){}
	virtual void reinit(const CalculatorParameterMap& params) = 0;
	virtual double eval(const CalculatorArgument * arg) = 0;
protected:
};

/* setups all values of cParameters having names 
 * from fParameters equal to their values
 */
void mergeParameters(CalculatorParameterMap& cParameters,
                     const FitParameterList& fParameters);

class NonlinearFitter
{
public:
	enum FitType {fitLIN, fitLOG};
	NonlinearFitter();
	virtual void init(FitCalculator * calculator,
			const FitParameterList& fitParameters,
			const CalculatorParameterMap& calcParameters,
			const DataPointList& dataPoints);
	virtual size_t getNbFitParameters() const {return m_FitParameters.size();}
	virtual const FitParameterList& getFitParameters() const {return m_FitParameters;}
	virtual const FitParameter& getFitParameter(size_t i) const {return m_FitParameters[i];}
	virtual void fit(FitType type = fitLIN, int nbit = 100) = 0;
	virtual ~NonlinearFitter();
	/*post fit information*/
	virtual ResidualValueType getFinit() const {return m_Finit;}
	virtual ResidualValueType getFfin() const {return m_Ffin;}
	virtual size_t getNbFuncEval() const {return m_nbFuncEval;}
	virtual size_t getNbGradEval() const {return m_nbGradEval;}
	virtual size_t getNbIterations() const {return m_nbIterPerf;}
	virtual std::string getName() const {return m_Name;}
	virtual double getCovarMatrix(size_t i, size_t j) const;
	virtual std::string getReasonToStop() const = 0;
	virtual int getReasonToStopId() const {return m_reasonToStopId;}
protected:
	/*resets CalculatorParameters according to x-array */
	virtual void resetCalculatorParameters(const double * x);
	/* resets FitParameters according to x-array */
	virtual void resetFitParameters(const double * x);

	/*map of calculator parameters*/
	CalculatorParameterMap m_CalculatorParameters;
	/*vector of fit parameters*/
	FitParameterList m_FitParameters;
	/*fit calculator*/
	FitCalculator * m_Calculator;
	/*'experimental' data points*/
	DataPointList m_DataPoints;

	/*post fit info*/
	ResidualValueType m_Finit, m_Ffin;
	size_t m_nbFuncEval;
	size_t m_nbGradEval;
	size_t m_nbIterPerf;
	int m_reasonToStopId;
	/*covariance matrix*/
	double ** m_matrixCovar;

	std::string m_Name;

	void allocateCovarMatrix();
	void deallocateCovarMatrix();
};
}
#endif /* NONLINEARFIT_H_ */
