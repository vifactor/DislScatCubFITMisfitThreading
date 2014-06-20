/*
 * PortFitter.cpp
 *
 *  Created on: 26 ���. 2013
 *      Author: kopp
 */

#include "PortFitter.h"

namespace NonlinearFit
{

PortFitter::PortFitter() : NonlinearFitter()
{
	m_Name = "PORT";
}

PortFitter::~PortFitter()
{
}

int LinNormFunc(int * n, int * m, double * x, int * nf, double * v,
		int * ui, double * ur, int (* uf)())
{
	register int i;
	static ResidualValueType calcResidual, expResidual;
	static PortFitter * fitter;

	/* ui is a pointer to the fitter */
	fitter = reinterpret_cast<PortFitter *>(ui);

	/* reinitialize calculator */
	fitter->resetCalculator(x);

	/* calculate residuals */
	for (i = 0; i < *n; ++i)
	{
		calcResidual = fitter->getCalcResidual(i);
		expResidual = fitter->getExpResidual(i);

		v[i] = calcResidual - expResidual;
	}

	return 0;
}

int LogNormFunc(int * n, int * m, double * x, int * nf, double * v,
		int * ui, double * ur, int (* uf)())
{
	register int i;
	static ResidualValueType calcResidual, expResidual;
	static PortFitter * fitter;

	/* ui is a pointer to the fitter */
	fitter = reinterpret_cast<PortFitter *>(ui);

	/* reinitialize calculator */
	fitter->resetCalculator(x);

	/* calculate residuals */
	for (i = 0; i < *n; ++i)
	{
		calcResidual = fitter->getCalcResidual(i);
		expResidual = fitter->getExpResidual(i);
		if(calcResidual > 0)
		{
			v[i] = log(calcResidual) - expResidual;
		}
		else
		{
			v[i] = 0. - expResidual;
		}
	}
	return 0;
}

void PortFitter::resetCalculator(const double * x)
{
	resetCalculatorParameters(x);
	m_Calculator->reinit(m_CalculatorParameters);
}

ResidualValueType PortFitter::getExpResidual(size_t i)
{
	return m_DataPoints[i].m_Residual;
}

ResidualValueType PortFitter::getCalcResidual(size_t i)
{
	return m_Calculator->eval(m_DataPoints[i].m_Argument);
}

void PortFitter::fit(FitType type, int nbit)
{

	//the value of f(xinit) for final output
	int nbParams, nbResids, liv, lv;
	double *x_ptr, *xb_ptr, *v_ptr;
	int *iv_ptr;
	int kind;
	int (* NormFunc)(int * m, int * n, double * x, int * nf, double * v,int * ui, double * ur, int (* uf)());

	nbParams = m_FitParameters.size();
	nbResids = m_DataPoints.size();

	if(type == fitLIN)
	{
		NormFunc = LinNormFunc;
	}else if(type == fitLOG)
	{
		NormFunc = LogNormFunc;
		for(size_t i = 0; i < m_DataPoints.size(); ++i)
		{
			if(m_DataPoints[i].m_Residual > 0)
			{
				m_DataPoints[i].m_Residual = log(m_DataPoints[i].m_Residual);
			}
			else
			{
				m_DataPoints[i].m_Residual = 0.0;
			}
		}
	}
	else
	{
		NormFunc = LinNormFunc;
	}

	// RN2GB: "LIV...... LENGTH OF IV... LIV MUST BE AT LEAST 4*P + 82"
	liv	= 4 * nbParams + 82;
	// RN2GB: "LV....... LENGTH OF V...  LV  MUST BE AT LEAST 105 + P*(N + 2*P + 17) + 2*N"
	lv = 105 + nbParams * (nbResids + 2 * nbParams + 17) + 2 * nbResids;

	x_ptr = new double[nbParams];
	xb_ptr = new double[2 * nbParams];
	v_ptr = new double[lv];
	iv_ptr = new int[liv];

	//initialize nl2sol parameters
	for (size_t i = 0; i < m_FitParameters.size(); ++i)
	{
		x_ptr[i] = m_FitParameters[i].m_Value;
		xb_ptr[2 * i] = m_FitParameters[i].m_Lbvalue;
		xb_ptr[2 * i + 1] = m_FitParameters[i].m_Ubvalue;
	}
	//initialize default values for iv and v optimization settings
	kind = 1;
	iv_ptr[0] = 0;
	divset_(&kind, iv_ptr, &liv, &lv, v_ptr);
	//turn of output
	iv_ptr[ivPRUNIT] = 0;
	//max iterations allowed
	iv_ptr[ivMXITER] = 0;
	//compute a covariance matrix !!!DN2FB DOES NOT COMPUTE THE COVARIANCE MATRIX!!!
	//iv_ptr.get()[ivRDREQ] = 1;

	//try to execute optimization algorithm
	dn2fb_(&nbResids, &nbParams, x_ptr, xb_ptr, NormFunc,
			iv_ptr, &liv, &lv, v_ptr,
			reinterpret_cast<int *>(this), NULL, NULL);

	//storage size for iv_ptr or v_ptr was too small
	if((iv_ptr[0] == 15) || (iv_ptr[0] == 16))
	{
		//reallocate iv_ptr and v_ptr
		liv = iv_ptr[43];
		lv = iv_ptr[44];
		delete[] iv_ptr;
		delete[] v_ptr;
		iv_ptr = new int[liv];
		v_ptr = new double[lv];
		//reset defaults
		divset_(&kind, iv_ptr, &liv, &lv, v_ptr);
		//turn of output
		iv_ptr[ivPRUNIT] = 0;
		//max iterations allowed
		iv_ptr[ivMXITER] = 0;

		//execute optimization algorithm once to get f(xinit)
		dn2fb_(&nbResids, &nbParams, x_ptr, xb_ptr, NormFunc,
				iv_ptr, &liv, &lv, v_ptr,
				reinterpret_cast<int *>(this), NULL, NULL);
		// f (x) at the start of the last iteration
		m_Finit = v_ptr[vF0];

		//max iterations allowed
		iv_ptr[ivMXITER] = nbit;
		//compute a covariance matrix
		//iv_ptr.get()[ivRDREQ] = 1;
		//execute optimization algorithm nbit-1 times to optimize
		dn2fb_(&nbResids, &nbParams, x_ptr, xb_ptr, NormFunc,
				iv_ptr, &liv, &lv, v_ptr,
				reinterpret_cast<int *>(this), NULL, NULL);
	}
	/*reset v(vF0) to m_Finit in order to have correct information at the end*/
	v_ptr[vF0] = m_Finit;
	resetFitParameters(x_ptr);
	extractFitInfo(iv_ptr, v_ptr);

	/*deallocate arrays*/
	delete[] iv_ptr;
	delete[] v_ptr;
	delete[] x_ptr;
	delete[] xb_ptr;
}

void PortFitter::extractFitInfo(const int* iv, const double* v)
{
	//size_t nbParams;

	//nbParams = m_FitParameters.size();
	m_Finit = sqrt(2 * v[vF0]);
	m_Ffin = sqrt(2 * v[vF]);
	m_nbFuncEval = iv[ivNFCALL] + iv[ivNGCALL];
	m_nbGradEval = iv[ivNGCALL] - iv[ivNGCOV];
	m_nbIterPerf = iv[ivNITER];
	m_reasonToStopId = iv[0];
	/*extract covariance matrix*/

	std::cout << "This routine does not compute the covariance matrix!" << std::endl;
}

std::string PortFitter::getReasonToStop() const
{
	switch(m_reasonToStopId)
	{
	case 3:
		return "X-convergence";
		break;
	case 4:
		return "Relative function convergence";
		break;
	case 5:
		return "Both X- and relative function convergence";
		break;
	case 6:
		return "Absolute function convergence";
		break;
	case 10:
		return "Max nb iterations";
		break;
	default:
		return "Invalid parameters";
		break;
	}
	return "Invalid parameters";
}

} /* namespace NonlinearFit */
