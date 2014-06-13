/*
 * MisfitInterfaceHex.cpp
 *
 *  Created on: 7 june 2013
 *      Author: kopp
 */

#include "MisfitInterfaceHex.h"

MisfitInterfaceHex::MisfitInterfaceHex(double rho, double bx, double bz, double Qx, double Qz,
		double nu, double d)
{
	alpha = 1.0 / (2 * (1- nu));
	alpha2 = alpha * alpha;

	bx2 = bx * bx;
	bz2 = bz * bz;

	Qx2 = Qx * Qx;
	Qz2 = Qz * Qz;
	QxQz = Qx * Qz;

	m_d = d;
	m_rho = rho;

	/*dummy initialization*/
	init(1.0);
}

MisfitInterfaceHex::~MisfitInterfaceHex()
{
}

double MisfitInterfaceHex::T(double x, double z1, double z2) const
{
	static double z, result;

	if (m_rho != 0.0)
	{
		init(z1);
		z = z1 - z2;

		/*
		 * division coefficient 3 appears due to 3 dislocation sets
		 * it is written here in order to equalize "real" Monte-Carlo dislocation density
		 *  and
		 * density present here
		*/
		result = m_rho / 2 / 3
				* (wxx() * x * x + wxz() * x * z + wzz() * z * z);
	}
	else
	{
		result = 0.0;
	}

	return result;
}


double MisfitInterfaceHex::wxx() const
{
	static double result, result_bx, result_bz;

	if (bx2 != 0.0)
	{
		result_bx = 3.0 / 8
				* (3 * Qx2 * Wxxxx_bx()
						+ 4 * Qz2 * Wzxzx_bx());
	}
	else
	{
		result_bx = 0.0;
	}

	if (bz2 != 0.0)
	{
		result_bz = 3.0 / 8
				* (3 * Qx2 * Wxxxx_bz()
						+ 4 * Qz2 * Wzxzx_bz());
	}
	else
	{
		result_bz = 0.0;
	}
	result = result_bx + result_bz;

	return result;
}

double MisfitInterfaceHex::wzz() const
{
	double result, result_bx, result_bz;

	if (bx2 != 0.0)
	{
		result_bx = 3.0 / 2 * Qx2 * Wxzxz_bx()
				+ 3 * Qz2 * Wzzzz_bx();
	}
	else
	{
		result_bx = 0.0;
	}

	if (bz2 != 0.0)
	{
		result_bz = 3.0 / 2 * Qx2 * Wxzxz_bz()
				+ 3 * Qz2 * Wzzzz_bz();
	}
	else
	{
		result_bz = 0.0;
	}

	result = result_bx + result_bz;
	return result;
}

double MisfitInterfaceHex::wxz() const
{
	double result, result_bx, result_bz;

	if (bx2 != 0.0)
	{
		result_bx = 3 * QxQz * (Wxxzz_bx() + Wxzzx_bx());
	}
	else
	{
		result_bx = 0.0;
	}

	if (bz2 != 0.0)
	{
		result_bz = 3 * QxQz * (Wxxzz_bz() + Wxzzx_bz());
	}
	else
	{
		result_bz = 0.0;
	}

	result = result_bx + result_bz;
	return result;
}

void MisfitInterfaceHex::init(double z) const
{
	xi = z / m_d;
	xi2 = xi * xi;
	xi3 = xi2 * xi;
	xi4 = xi3 * xi;
	xi5 = xi4 * xi;
	xi6 = xi5 * xi;
	xi7 = xi6 * xi;
	xi8 = xi7 * xi;

	denom = (8. * m_d * M_PI * (-1 + xi) * gsl_pow_5(1 + xi));
}

double MisfitInterfaceHex::Wxxxx_bx() const
{
	static double result;

	result = bx2
			* (-2 - 6 * xi - 2 * (5 + 2 * alpha) * xi2
					+ 2 * (-5 + 4 * alpha) * xi3
					+ (-4 + 18 * alpha - 21 * alpha2) * xi4
					+ (8 - 9 * alpha) * alpha * xi5
					+ alpha * (2 + 5 * alpha) * xi6 + 7 * alpha2 * xi7
					+ 2 * alpha2 * xi8) / denom;

	//std::cout << "Wxxxx_bx:\t" << result << std::endl;
	return result;
}

double MisfitInterfaceHex::Wzzzz_bx() const
{
	double result;

	result = bx2
			* (-2 * gsl_pow_2(1 + xi) * (1 + xi + 2 * xi2)
					- 2 * alpha
							* (-4 - 12 * xi - 22 * xi2 - 16 * xi3 + xi4
									+ 4 * xi5 + xi6)
					+ alpha2
							* (-8 - 24 * xi - 48 * xi2 - 24 * xi3 - xi4
									+ 7 * xi5 + 9 * xi6 + 7 * xi7 + 2 * xi8))
			/ denom;

	//std :: cout << "Wzzzz_bx:\t" << result << std::endl;
	return result;
}

double MisfitInterfaceHex::Wxzxz_bx() const
{
	double result;

	result = -bx2
			* (2 + (6 - 4 * alpha) * xi
					+ 2 * (5 + 2 * alpha + 4 * alpha2) * xi2
					+ (10 + 12 * alpha - 4 * alpha2) * xi3
					+ (4 + 10 * alpha - 3 * alpha2) * xi4
					- (-8 + alpha) * alpha * xi5
					+ alpha * (2 + 7 * alpha) * xi6 + 7 * alpha2 * xi7
					+ 2 * alpha2 * xi8) / denom;
	//std :: cout << "Wxzxz_bx:\t" << result << std::endl;
	return result;
}


double MisfitInterfaceHex::Wxxzz_bx() const
{
	double result;

	result = -bx2
			* (2 * gsl_pow_2(1 + xi) * (1 + xi + 2 * xi2)
					- 4 * alpha * gsl_pow_2(1 + xi) * (1 + xi + 2 * xi2)
					+ alpha2 * xi2
							* (-4 + 8 * xi - 3 * xi2 - xi3 + 7 * xi4 + 7 * xi5
									+ 2 * xi6)) / denom;
	//std :: cout << "Wxxzz_bx:\t" << result << std::endl;
	return result;
}

double MisfitInterfaceHex::Wzxzx_bx() const
{
	double result;

	result = -bx2
			* (2 + (6 + 4 * alpha) * xi
					+ 2 * (5 - 2 * alpha + 4 * alpha2) * xi2
					- 2 * (-5 + 6 * alpha + 2 * alpha2) * xi3
					+ (4 - 10 * alpha - 3 * alpha2) * xi4
					- alpha * (8 + alpha) * xi5
					+ alpha * (-2 + 7 * alpha) * xi6 + 7 * alpha2 * xi7
					+ 2 * alpha2 * xi8) / denom;

	//std :: cout << "Wzxzx_bx:\t" << result << std::endl;
	return result;
}

double MisfitInterfaceHex::Wxzzx_bx() const
{
	double result;

	result = bx2
			* (2 + 6 * xi + (10 - 8 * alpha2) * xi2 + 2 * (5 + 2 * alpha2) * xi3
					+ (4 + 3 * alpha2) * xi4 + alpha2 * xi5 - 7 * alpha2 * xi6
					- 7 * alpha2 * xi7 - 2 * alpha2 * xi8) /denom;

	//std :: cout << "Wxzzx_bx:\t" << result << std::endl;
	return result;
}

double MisfitInterfaceHex::Wxxxx_bz() const
{

	double result;

	result = bz2
			* (-2 + 2 * xi - 2 * (-5 + 6 * alpha) * xi2 + (-2 + 8 * alpha) * xi3
					- 5 * (4 - 10 * alpha + 9 * alpha2) * xi4
					+ (-16 + 40 * alpha - 45 * alpha2) * xi5
					+ (-4 + 10 * alpha - 5 * alpha2) * xi6 + 11 * alpha2 * xi7
					+ 4 * alpha2 * xi8) / denom;
	//std :: cout << "Wxxxx_bz:\t" << result << std::endl;
	return result;
}

double MisfitInterfaceHex::Wzzzz_bz() const
{
	double result;

	result = bz2
			* (-2 * gsl_pow_2(1 + xi) * (1 - 3 * xi + 4 * xi3 + 2 * xi4)
					+ alpha
							* (8 - 8 * xi - 28 * xi2 + 30 * xi4 + 24 * xi5
									+ 6 * xi6)
					+ alpha2
							* (-8 + 8 * xi + 16 * xi2 + 8 * xi3 - 25 * xi4
									- 29 * xi5 - xi6 + 11 * xi7 + 4 * xi8))
			/ denom;
	//std :: cout << "Wzzzz_bz:\t" << result << std::endl;
	return result;
}

double MisfitInterfaceHex::Wxzxz_bz() const
{
	double result;

	result = -bz2
			* (10 + (30 - 20 * alpha) * xi
					+ (30 - 28 * alpha + 24 * alpha2) * xi2
					+ 2 * (1 - 2 * alpha + 6 * alpha2) * xi3
					- 5 * (4 - 2 * alpha + 3 * alpha2) * xi4
					+ (-16 + 8 * alpha - 21 * alpha2) * xi5
					+ (-4 + 2 * alpha + alpha2) * xi6 + 11 * alpha2 * xi7
					+ 4 * alpha2 * xi8) / denom;
	//std :: cout << "Wxzxz_bz:\t" << result << std::endl;
	return result;
}

double MisfitInterfaceHex::Wxxzz_bz() const
{
	double result;

	result = -bz2
			* (2 * gsl_pow_2(1 + xi) * (1 - 3 * xi + 4 * xi3 + 2 * xi4)
					- 4 * alpha * gsl_pow_2(1 + xi)
							* (1 - 3 * xi + 4 * xi3 + 2 * xi4)
					+ alpha2 * xi2
							* (-12 + 8 * xi + 5 * xi2 - 5 * xi3 + 5 * xi4
									+ 11 * xi5 + 4 * xi6)) / denom;
	//std :: cout << "Wxxzz_bz:\t" << result << std::endl;
	return result;
}

double MisfitInterfaceHex::Wzxzx_bz() const
{
	double result;

	result = -bz2
			* (10 + 10 * (3 + 2 * alpha) * xi
					+ (30 + 28 * alpha + 24 * alpha2) * xi2
					+ 2 * (1 + 2 * alpha + 6 * alpha2) * xi3
					- 5 * (4 + 2 * alpha + 3 * alpha2) * xi4
					- (16 + 8 * alpha + 21 * alpha2) * xi5
					+ (-4 - 2 * alpha + alpha2) * xi6 + 11 * alpha2 * xi7
					+ 4 * alpha2 * xi8) / denom;
	//std :: cout << "Wzxzx_bz:\t" << "("<< xi <<"):\t"  << result << std::endl;
	return result;
}

double MisfitInterfaceHex::Wxzzx_bz() const
{
	double result;

	result = -(bz2
			* (-10 - 30 * xi + 6 * (-5 + 4 * alpha2) * xi2
					+ 2 * (-1 + 6 * alpha2) * xi3 - 5 * (-4 + 3 * alpha2) * xi4
					+ (16 - 21 * alpha2) * xi5 + (4 + alpha2) * xi6
					+ 11 * alpha2 * xi7 + 4 * alpha2 * xi8)) / denom;
	//std :: cout << "Wxzzx_bz:\t" << result << std::endl;

	return result;
}
