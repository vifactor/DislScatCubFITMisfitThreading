/*
 * ANASampleCub.h
 *
 *  Created on: 12 черв. 2013
 *      Author: kopp
 */

#ifndef ANASAMPLECUB_H_
#define ANASAMPLECUB_H_

#include "ANAMisfitInterfaceCub.h"
#include "ANAThreadingLayerCub.h"
#include <vector>

class ANASampleCub
{
public:
	ANASampleCub(double thickness, double size);
	~ANASampleCub();

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
	std::vector<ANAMisfitInterfaceCub * > m_interfaces;
	std::vector<ANAThreadingLayerCub * > m_layers;
};

#endif /* ANASampleCub_H_ */
