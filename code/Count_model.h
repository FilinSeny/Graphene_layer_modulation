#ifndef COUNT_MODEL_H_INCLUDED
#define COUNT_MODEL_H_INCLUDED




#endif // COUNT_MODEL_H_INCLUDED

#pragma once

#include<iostream>
#include<cmath>
#include<fstream>
#include <complex>
#include <ostream>
#include <math.h>

const long double PI = 3.14159265358979323846;
const long double V_Light = 299792458.;
const long double eps0_vac = 8.85418781762039 * 1e-12;
const long double mu0_vac = 4. * PI * 1e-7;
const std::complex<long double> i(0, 1);

const std::complex<long double> r1(1., 0.); ///для вакуума очев
const std::complex<long double> r2(1.4504, 0.); ///SiO2 при 1мкм
const std::complex<long double> r3(1.436, 9.495); ///AL при 1 мкм

const std::complex<long double> eps1 = r1 * r1 * eps0_vac;
const std::complex<long double> eps2 = r2 * r2 * eps0_vac;
const std::complex<long double> eps3 = r3 * r3 * eps0_vac;

const long double vawe_len = 1e-6; //длинна волны в вакууме
const long double omega = 2. * PI / vawe_len / sqrt(eps0_vac * mu0_vac);

const std::complex<long double> k_1 = omega * sqrt(eps0_vac * mu0_vac);
const std::complex<long double> k_2 = k_1 * r2;
const std::complex<long double> k_3 = k_1 * r3;

const std::complex<long double> k_12 = k_1 * k_1;
const std::complex<long double> k_22 = k_2 * k_2;
const std::complex<long double> k_32 = k_3 * k_3;

const std::complex<long double> sigG = std::complex<long double>(1.0E-8, 1.0E-5) / eps0_vac;






class Count_Model {
public:
	const long double z_0;


	Count_Model(long double pos);

	std::complex<long double> Count_Hankel_Tranform(std::complex<long double>(*func)(const Count_Model&, const long double &, const long double &, const long double &), ///преобразуемая функция,
		const long double & z, const long double &R, double Bessel_num) const;

	friend std::complex<long double> U1(const Count_Model&, const long double & ,const long double &,const long double &);
	friend std::complex<long double> U0(const Count_Model&, long double, const long double, const long double);
	friend std::complex<long double> U(const Count_Model&, const long double &, const long double &, const long double &);
	friend std::complex<long double> W(const Count_Model&, const long double &, const long double &, const long double &);
	friend std::complex<long double> WP(const Count_Model&, const long double &, const long double &, const long double &);
	friend std::complex<long double> WP1(const Count_Model&, const long double &, const long double &, const long double &);
	friend std::complex<long double> WP_1(const Count_Model&, const long double &, const long double &, const long double &);

	friend std::complex<long double> d2zrF(const Count_Model&, const long double & rho, const long double & z, const long double & z0);
	friend std::complex<long double> d2zzF(const Count_Model&, const long double & rho, const long double & z, const long double & z0);
	friend std::complex<long double> dzF(const Count_Model&, const long double & rho, const long double & z, const long double & z0);
	friend std::complex<long double> drF(const Count_Model&, const long double & rho, const long double & z, const long double & z0);
	friend std::complex<long double> F(const Count_Model&, const long double & rho, const long double & z, const long double & z0);

	friend std::complex<long double> E_rho(const Count_Model&, const long double &, const long double &, const long double &);
	friend std::complex<long double> E_phi(const Count_Model&, const long double &, const long double &, const long double &);
	friend std::complex<long double> E_z(const Count_Model&, const long double &, const long double &, const long double &);

	friend std::complex<long double> H_rho(const Count_Model&, const long double &, const long double &, const long double &);
	friend std::complex<long double> H_phi(const Count_Model&, const long double &, const long double &, const long double &);
	friend std::complex<long double> H_z(const Count_Model&, const long double &, const long double &, const long double &);

	friend std::complex<long double> F_wg(const Count_Model&, const long double &, const long double &);
    friend std::complex<long double> F_ug(const Count_Model&, const long double &, const long double &);

	friend std::ostream & operator << (std::ostream & out, const Count_Model&);



};
