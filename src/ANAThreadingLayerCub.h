/*
 * ANAThreadingLayerCub.h
 *
 *  Created on: 11 черв. 2013
 *      Author: kopp
 */

#ifndef ANATHREADINGLAYERCUB_H_
#define ANATHREADINGLAYERCUB_H_

#include <gsl/gsl_math.h>
#include <iostream>

class ANAThreadingLayerCub
{
public:
	ANAThreadingLayerCub(double rho, double b_edge, double b_screw, double rc, double Qx, double Qz, double nu);
	virtual ~ANAThreadingLayerCub();

	double T(double r, double phi) const;
	double m_rho;
	double m_Rc;
protected:
	double m_gb_screw, m_gb2_screw;
	double m_gb_edge, m_gb2_edge;

	double m_C_screw, m_C1_edge, m_C2_edge;

	double T_screw(double r) const;
	double T_edge(double r, double phi) const;
	inline double chi_screw() const;
	inline double chi_edge(double phi) const;
};

#endif /* ANAThreadingLayerCub_H_ */
