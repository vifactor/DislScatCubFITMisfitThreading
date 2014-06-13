/*
 * ThreadingLayerHex.cpp
 *
 *  Created on: 11 черв. 2013
 *      Author: kopp
 */

#include "ThreadingLayerHex.h"

ThreadingLayerHex::ThreadingLayerHex(double rho, double b_edge, double b_screw, double rc, double Qx, double Qz, double nu)
{
	double Q, sin2Psi, cos2Psi;

	Q = sqrt(Qx * Qx + Qz * Qz);
	sin2Psi = gsl_pow_2(Qz / Q);
	cos2Psi = 1 - sin2Psi;

	m_C_screw = M_PI * sin2Psi;
	m_C1_edge = M_PI * (9 - 16 * nu + 8 * nu * nu) * cos2Psi
			/ (8 * gsl_pow_2(1 - nu));
	m_C2_edge = M_PI * (-2 * (3 - 4 * nu)) * cos2Psi
			/ (8 * gsl_pow_2(1 - nu));

	m_rho = rho;
	m_Rc = rc;

	m_gb2_edge = gsl_pow_2(Q * b_edge / (2 * M_PI));
	m_gb2_screw = gsl_pow_2(Q * b_screw / (2 * M_PI));

	m_gb_edge = sqrt(m_gb2_edge);
	m_gb_screw = sqrt(m_gb2_screw);
}

ThreadingLayerHex::~ThreadingLayerHex()
{
}

double ThreadingLayerHex::T(double r, double phi) const
{
	static double result;

	result = T_edge(r, phi) + T_screw(r);

	return result;
}

double ThreadingLayerHex::T_edge(double r, double phi) const
{
	static double result;

	if((m_gb2_edge > 0) && (m_rho > 0))
		result = 0.5 * m_rho * m_gb2_edge * r * r * chi_edge(phi)
				* log((m_Rc / m_gb_edge + r) / r);
	else
		result = 0.0;

	return result;
}

double ThreadingLayerHex::T_screw(double r) const
{
	static double result;

	if((m_gb2_screw > 0)&&(m_rho > 0))
		result = 0.5 * m_rho * m_gb2_screw * r * r * chi_screw()
				* log((m_Rc / m_gb_screw + r) / r);
	else
		result = 0.0;

	return result;
}

double ThreadingLayerHex::chi_screw() const
{
	return m_C_screw;
}

double ThreadingLayerHex::chi_edge(double phi) const
{
	return m_C1_edge + m_C2_edge * gsl_pow_2(cos(phi));
}
