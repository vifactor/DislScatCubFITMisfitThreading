/*
 * ANACalculatorCoplanarTriple.cpp
 *
 *  Created on: 31 june 2013
 *      Author: kopp
 */

#include "ANACalculatorCoplanarTriple.h"

double ana_coplanar_triple_integrand_xz1(double x, void *params)
{
	ANACalculatorCoplanarTriple * calculator;
	static double result;

	calculator = static_cast<ANACalculatorCoplanarTriple *> (params);

	result = exp(-calculator->m_coefficient * x * x - calculator->T_threading(x))
			* cos(calculator->m_frequency * x);

	return result;
}

ANACalculatorCoplanarTriple::ANACalculatorCoplanarTriple(size_t sampling, double epsabs)
{

	m_sample = NULL;

	m_qx = 0.0;
	m_qz = 0.0;

	m_resol2_x = 0.0;
	m_resol2_z = 0.0;

	m_scale = 0.0;
	m_coefficient = 0.0;
	m_frequency = 0.0;

	/*initialization of integration staff*/
	m_limit = 10000;
	m_epsabs = epsabs;
	m_function.function = &ana_coplanar_triple_integrand_xz1;
	m_function.params = this;

	m_qawo_table = gsl_integration_qawo_table_alloc(0.0, 0.0, GSL_INTEG_COSINE,
			m_limit);
	m_workspace = gsl_integration_workspace_alloc(m_limit);
	m_cyclic_workspace = gsl_integration_workspace_alloc(m_limit);

	m_a = 0.0;

	/*allocate and initialize z-arrays*/
	m_sampling = sampling;
	m_z1 = new double[m_sampling];
	m_integrand_values = new double[m_sampling];

	/*allocate interpolation staff*/
	m_interp =  gsl_interp_alloc (gsl_interp_cspline, m_sampling);
	m_accel = gsl_interp_accel_alloc ();

	gsl_set_error_handler_off ();
}

ANACalculatorCoplanarTriple::~ANACalculatorCoplanarTriple()
{
	delete[] m_z1;
	delete[] m_integrand_values;

	gsl_interp_free (m_interp);
	gsl_interp_accel_free(m_accel);

	if(m_qawo_table)
		gsl_integration_qawo_table_free(m_qawo_table);
	if(m_workspace)
		gsl_integration_workspace_free(m_workspace);
	if(m_cyclic_workspace)
		gsl_integration_workspace_free(m_cyclic_workspace);
}

void ANACalculatorCoplanarTriple::setSample(ANASampleCub * sample)
{
    m_sample = sample;
    for(size_t i = 0; i < m_sampling - 1; ++i)
    {
        m_z1[i] = i * m_sample->m_thickness / m_sampling;
        m_integrand_values[i] = 0.0;
    }
    /*at z1 = d, the correlation function is undefined*/
    m_z1[m_sampling - 1] = m_sample->m_thickness * 0.999;
}

void ANACalculatorCoplanarTriple::setResolution(double fwhm_qx, double fwhm_qz)
{
	/* FWHM = 2 * sqrt(2 * log(2)) / sigma */
	/* resol2 = 1 / (2 * sigma**2) */
	m_resol2_x = gsl_pow_2(fwhm_qx / 4) / log(2.0);
	m_resol2_z = gsl_pow_2(fwhm_qz / 4) / log(2.0);
}

double ANACalculatorCoplanarTriple::T_threading(double x) const
{
	static double result;

	result = m_sample->T_threading(fabs(x), 0.0);

	return result;
}

double
ANACalculatorCoplanarTriple::I(const double qx, const double qz, double & err) const
{
	static double result, abserr;

	m_qx = qx;
	m_qz = qz;

	err = 0.0;
	/*here integration over x is performed*/
	for(size_t i = 0; i < m_sampling; ++i)
	{
		result = 0.0;
		abserr = 0.0;
		prepare(m_z1[i]);
		//gsl_integration_qawo_table_set(m_qawo_table, m_frequency, 1, GSL_INTEG_COSINE);
		//gsl_integration_qawf (&m_function, m_a, m_epsabs, m_limit, m_workspace, m_cyclic_workspace, m_qawo_table, &result, &abserr);

		gsl_integration_qagiu (&m_function, m_a, m_epsabs, 1.e-6, m_limit, m_workspace, &result, &abserr);
		m_integrand_values[i] = 2 * result * m_scale;
		err += 2 * abserr;

		//std::cout << m_z1[i] << "\t" << m_integrand_values[i] << std::endl;
	}

	/*here integration over z1 is performed*/
	//double gsl_interp_eval_integ (const gsl_interp * interp, const double xa[], const double ya[], double a, double b, gsl_interp_accel * acc)
	gsl_interp_init (m_interp, m_z1, m_integrand_values, m_sampling);
	result = gsl_interp_eval_integ (m_interp, m_z1, m_integrand_values, m_z1[0], m_z1[m_sampling-1], m_accel);

	return result;
}

double ANACalculatorCoplanarTriple::I(const double qx, const double qz) const
{
	static double result, abserr;

	result = I(qx, qz, abserr);

	return result;
}

void ANACalculatorCoplanarTriple::prepare(double z1) const
{
	static double wxx, wxz, wzz;

	m_sample->wij(z1, wxx, wxz, wzz);

	/*take into account resolution*/
	wxx += m_resol2_x;
	wzz += m_resol2_z;

	m_scale = sqrt(M_PI / wzz) * exp(-m_qz * m_qz / (4 * wzz));
	m_coefficient = wxx - wxz * wxz / (4 * wzz);
	m_frequency = m_qx - m_qz * wxz / (2 * wzz);
}
