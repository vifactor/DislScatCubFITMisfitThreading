/*
 * PortFitter.h
 *
 *  Created on: 26 ���. 2013
 *      Author: kopp
 */

#ifndef PORTFITTER_H_
#define PORTFITTER_H_

#include "NonlinearFit.h"
#include <port.h>

namespace NonlinearFit
{

class PortFitter: public NonlinearFit::NonlinearFitter
{
public:
	PortFitter();
	virtual void fit(FitType type = fitLIN, int nbit=100);
	virtual std::string getReasonToStop() const;
	virtual ~PortFitter();
protected:
	enum IVIndices
	{
		ivNFCALL = 5, /*iv(ivNFCALL) nb of function calculations performed*/
		ivCOVPRT = 13,/*controls printing of a covariance matrix and regression diagnostic array*/
		ivNITER = 30, /*iv(ivNITER) nb of iterations performed*/
		ivMXITER = 17, /*iv(ivMXITER) the maximum number of iterations allowed*/
		ivPRUNIT = 20, /*iv(ivPRUNIT) the Fortran output unit for printing*/
		ivCOVMAT = 25, /*iv(ivCOVMAT) if positive, is the starting subscript in V for the lower triangle of covariance matrix*/
		ivNGCALL = 29, /*iv(ivNFCALL) is either the total number of gradient evaluations performed (when you provide
							analytic derivatives) or the number of additional function evaluations required to compute
							finite-difference derivative approximations*/
		ivNGCOV = 52, /* iv(ivNGCOV) is the number of additional gradient evaluations performed just for computing a covariance matrix or regression diagnostics*/
		ivRDREQ = 56 /*iv(ivRDREQ) tells whether to compute a covariance matrix or regression diagnostic array*/
	};
	enum VIndices
	{
		vF = 9, /*is the current function value f(x_fin)*/
		vF0 = 12/*v(vF0) is the function value of f (x) at the start of the last iteration*/
	};
	friend int LinNormFunc(int * n, int * m, double * x, int * nf, double * v,
			int * ui, double * ur, int (* uf)());
	friend int LogNormFunc(int * n, int * m, double * x, int * nf, double * v,
			int * ui, double * ur, int (* uf)());

	void resetCalculator(const double * x);
	ResidualValueType getExpResidual(size_t i);
	ResidualValueType getCalcResidual(size_t i);

	void extractFitInfo(const int* iv, const double* v);
};

} /* namespace NonlinearFit */
#endif /* PORTFITTER_H_ */
