/*
 * ANASampleCub.cpp
 *
 *  Created on: 12 june 2013
 *      Author: kopp
 */

#include "ANASampleCub.h"

ANASampleCub::ANASampleCub(double thickness, double size)
{
	m_thickness = thickness; m_size = size;
}

ANASampleCub::~ANASampleCub()
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

double ANASampleCub::T_threading(double r, double phi) const
{
	static double result;

	result = 0.0;
	for(size_t i = 0; i < m_layers.size(); ++i)
	{
		result += m_layers[i]->T(r, phi);
	}
	return result;
}

double ANASampleCub::T_misfit(double x, double y, double z1, double z2) const
{
	static double result;

	result = 0.0;
	for(size_t i = 0; i < m_interfaces.size(); ++i)
	{
		result += m_interfaces[i]->T(x, y, z1, z2);
	}
	return result;
}

void ANASampleCub::addThreadingLayer(double rho, double b_edge, double b_screw, double rc,
		double Qx, double Qz, double nu)
{
    m_layers.push_back(new ANAThreadingLayerCub(rho, b_edge, b_screw, rc, Qx, Qz, nu));
}

void ANASampleCub::addMisfitInterface(double rho, double bx, double by, double bz,
                    double Qx, double Qy, double Qz,
                    double phi, double nu, double d)
{
    m_interfaces.push_back(new AnalyticalMisfitInterfaceCub(rho, bx, by, bz,
                         Qx, Qy, Qz,
                         phi, nu, d));
}

void ANASampleCub::wij(double z1, double& wxx, double& wxz, double& wzz) const
{
	wxx = 0.0;
	wxz = 0.0;
	wzz = 0.0;

	for(size_t i = 0; i < m_interfaces.size(); ++i)
	{
		m_interfaces[i]->init(z1);

		wxx += m_interfaces[i]->wxx() * m_interfaces[i]->get_rho() / 2;
		wxz += m_interfaces[i]->wxz() * m_interfaces[i]->get_rho() / 2;
		wzz += m_interfaces[i]->wzz() * m_interfaces[i]->get_rho() / 2;
	}
	
	/*std::cout << "wxx:\t" << wxx << std::endl;
	std::cout << "wxz:\t" << wxz << std::endl;
	std::cout << "wzz:\t" << wzz << std::endl;*/
}

void ANASampleCub::resetMisfitInterface(size_t i, double rho)
{
	if(i < m_interfaces.size())
	{
		m_interfaces[i]->set_rho(rho);
	}
}

void ANASampleCub::resetThreadingLayer(size_t i, double rho, double rc)
{
	if(i < m_layers.size())
	{
		m_layers[i]->m_rho = rho;
		m_layers[i]->m_Rc = rc;
	}
}
