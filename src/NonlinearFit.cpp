/*
 * NonlinearFitter.cpp
 *
 *  Created on: 13 july 2012
 *      Author: kopp
 */

#include "NonlinearFit.h"

using namespace NonlinearFit;

NonlinearFitter::NonlinearFitter()
{
	m_Calculator = NULL;
	m_Finit = 0.0;
	m_Ffin = 0.0;
	m_nbFuncEval = 0;
	m_nbGradEval = 0;
	m_nbIterPerf = 0;
	m_reasonToStopId = 0;
	m_matrixCovar = NULL;
}

NonlinearFitter::~NonlinearFitter()
{
	if (m_matrixCovar != NULL)
	{
		deallocateCovarMatrix();
	}
}

void NonlinearFitter::resetCalculatorParameters(const double * x)
{
	for (size_t i = 0; i < m_FitParameters.size(); ++i)
	{
		m_CalculatorParameters[m_FitParameters[i].m_Name] = x[i];
	}
}

void NonlinearFitter::resetFitParameters(const double * x)
{
	for (size_t i = 0; i < m_FitParameters.size(); ++i)
	{
		m_FitParameters[i].m_Value = x[i];
	}
}

void NonlinearFitter::init(FitCalculator * calculator,
		const FitParameterList& fitParameters,
		const CalculatorParameterMap& calcParameters,
		const DataPointList& dataPoints)
{
	//size_t nbParams;

	m_Calculator = calculator;
	m_FitParameters = fitParameters;
	m_CalculatorParameters = calcParameters;
	m_DataPoints = dataPoints;

	//nbParams = m_FitParameters.size();

	if(m_matrixCovar != NULL)
	{
		deallocateCovarMatrix();
	}
	allocateCovarMatrix();
}

void NonlinearFit::mergeParameters(CalculatorParameterMap& cParameters,
		const FitParameterList& fParameters)
{
	for (size_t i = 0; i < fParameters.size(); ++i)
	{
		cParameters[fParameters[i].m_Name] = fParameters[i].m_Value;
	}
}

double NonlinearFitter::getCovarMatrix(size_t i, size_t j) const
{
	if(m_matrixCovar)
	{
		return m_matrixCovar[i][j];
	}
	else
	{
		return 0.0;
	}
}

void NonlinearFitter::allocateCovarMatrix()
{
	m_matrixCovar = new double *[m_FitParameters.size()];
	for (size_t i = 0; i < m_FitParameters.size(); ++i)
	{
		m_matrixCovar[i] = new double[m_FitParameters.size()];
	}
}

void NonlinearFitter::deallocateCovarMatrix()
{
	for (size_t i = 0; i < m_FitParameters.size(); ++i)
	{
		delete[] m_matrixCovar[i];
	}
	delete[] m_matrixCovar;
}
