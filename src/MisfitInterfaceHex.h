/*
 * TMisfitInterfaceHex.h
 *
 *  Created on: 7 черв. 2013
 *      Author: kopp
 */

#ifndef MISFITHEX_H_
#define MISFITHEX_H_

#include <gsl/gsl_math.h>
#include <iostream>

class MisfitInterfaceHex
{
public:
	MisfitInterfaceHex(double rho, double bx, double bz, double Qx, double Qz, double nu,
		 double d);
	virtual ~MisfitInterfaceHex();
    double wxx() const;
    double wzz() const;
    double wxz() const;
    void init(double z)const;
    double T(double x, double z1, double z2) const;

    double m_rho;
    double m_d;
protected:

    double Wxxxx_bx() const;
    double Wzzzz_bx() const;
    double Wxzxz_bx() const;
    double Wxxzz_bx() const;
    double Wzxzx_bx() const;
    double Wxzzx_bx() const;

    double Wxxxx_bz() const;
    double Wzzzz_bz() const;
    double Wxzxz_bz() const;
    double Wxxzz_bz() const;
    double Wzxzx_bz() const;
    double Wxzzx_bz() const;

    double alpha, alpha2;
    double bx2, bz2;
    double Qx2, Qz2, QxQz;
    mutable double xi, xi2, xi3, xi4, xi5, xi6, xi7, xi8;
    mutable double denom;
};

#endif /* MISFITHEX_H_ */
