/*
 * ANACalculatorCoplanarTriple.h
 *
 *  Created on: 21 june 2013
 *      Author: kopp
 */

#ifndef ANACALCULATORCOPLANARTRIPLE_H_
#define ANACALCULATORCOPLANARTRIPLE_H_

#include "ANASampleCub.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_errno.h>

class ANACalculatorCoplanarTriple
{
public:
	ANACalculatorCoplanarTriple(size_t sampling = 100, double epsabs = 1e-5);
	virtual ~ANACalculatorCoplanarTriple();

	virtual void setResolution(double fwhm_qx, double fwhm_qz);
	virtual void setSample(ANASampleCub * );
	ANASampleCub * getSample() {return m_sample;}
	virtual double I(const double qx, const double qz) const;
	virtual double I(const double qx, const double qz, double & error) const;
	
	friend double ana_coplanar_triple_integrand_xz1(double x, void *params);
protected:
    ANASampleCub * m_sample;

	mutable double m_qx, m_qz;
	double m_resol2_x, m_resol2_z;

	size_t m_sampling;
	mutable double * m_z1, * m_integrand_values;

	void prepare(double z1) const;
	inline double T_threading(double x) const;

	mutable double m_scale, m_coefficient, m_frequency;

	/*gsl integration staff*/
	gsl_integration_workspace * m_workspace;
	gsl_integration_workspace * m_cyclic_workspace;
	gsl_integration_qawo_table * m_qawo_table;
	mutable gsl_function m_function;
	size_t m_limit;
	double m_a, m_epsabs;

	/*gsl interpolation (and integration) staff*/
	gsl_interp * m_interp;
	gsl_interp_accel * m_accel;
};

#endif /* ANACALCULATORCOPLANARTRIPLE_H_ */
