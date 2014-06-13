/*
 * ANASampleHex.cpp
 *
 *  Created on: 12 june 2013
 *      Author: kopp
 */

#include "ANASampleHex.h"

ANASampleHex::ANASampleHex(double thickness, double size)
{
	m_thickness = thickness; m_size = size;
}

ANASampleHex::~ANASampleHex()
{
	for(size_t i = 0; i < m_interfaces.size(); ++i)
	{
		delete m_interfaces[i];
	}
	for(size_t i = 0; i < m_layers.size(); ++i)
	{
		delete m_layers[i];
	}
}

double ANASampleHex::T_threading(double r, double phi) const
{
	static double result;

	result = 0.0;
	for(size_t i = 0; i < m_layers.size(); ++i)
	{
		result += m_layers[i]->T(r, phi);
	}
	return result;
}

double ANASampleHex::T_misfit(double x, double z1, double z2) const
{
	static double result;

	result = 0.0;
	for(size_t i = 0; i < m_interfaces.size(); ++i)
	{
		result += m_interfaces[i]->T(x, z1, z2);
	}
	return result;
}

void ANASampleHex::addThreadingLayer(double rho, double b_edge, double b_screw, double rc,
		double Qx, double Qz, double nu)
{
	m_layers.push_back(new ThreadingLayerHex(rho, b_edge, b_screw, rc, Qx, Qz, nu));
}

void ANASampleHex::addMisfitInterface(double rho, double bx, double bz, double Qx,
		double Qz, double nu, double d)
{
	m_interfaces.push_back(new MisfitInterfaceHex(rho, bx, bz, Qx, Qz, nu, d));
}

void ANASampleHex::w_matrix(double z1, double& wxx, double& wxz, double& wzz) const
{
	wxx = 0.0;
	wxz = 0.0;
	wzz = 0.0;

	for(size_t i = 0; i < m_interfaces.size(); ++i)
	{
		m_interfaces[i]->init(z1);

		wxx += m_interfaces[i]->wxx() * m_interfaces[i]->m_rho / 2;
		wxz += m_interfaces[i]->wxz() * m_interfaces[i]->m_rho / 2;
		wzz += m_interfaces[i]->wzz() * m_interfaces[i]->m_rho / 2;
	}
}

void ANASampleHex::resetMisfitInterface(size_t i, double rho)
{
	if(i < m_interfaces.size())
	{
		m_interfaces[i]->m_rho = rho;
	}
}

void ANASampleHex::resetThreadingLayer(size_t i, double rho, double rc)
{
	if(i < m_layers.size())
	{
		m_layers[i]->m_rho = rho;
		m_layers[i]->m_Rc = rc;
	}
}
