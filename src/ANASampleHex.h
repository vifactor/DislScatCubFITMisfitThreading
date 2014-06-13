/*
 * ANASampleHex.h
 *
 *  Created on: 12 черв. 2013
 *      Author: kopp
 */

#ifndef ANASAMPLEHEX_H_
#define ANASAMPLEHEX_H_

#include "MisfitInterfaceHex.h"
#include "ThreadingLayerHex.h"
#include <vector>

class ANASampleHex
{
public:
	ANASampleHex(double thickness, double size);
	~ANASampleHex();

	void addThreadingLayer(double rho, double b_edge, double b_screw, double rc,
			double Qx, double Qz, double nu);
	void addMisfitInterface(double rho, double bx, double bz, double Qx,
			double Qz, double nu, double d);
	void resetMisfitInterface(size_t i, double rho);
	void resetThreadingLayer(size_t i, double rho, double rc);
	double T_threading(double r, double phi) const;
	double T_misfit(double x, double z1, double z2) const;
	void w_matrix(double z1, double& wxx, double& wxz, double& wzz) const;
	double m_thickness;
	double m_size;
protected:
	std::vector<MisfitInterfaceHex * > m_interfaces;
	std::vector<ThreadingLayerHex * > m_layers;
};

#endif /* ANASAMPLEHEX_H_ */
