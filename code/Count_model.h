#ifndef COUNT_MODEL_H_INCLUDED
#define COUNT_MODEL_H_INCLUDED




#endif // COUNT_MODEL_H_INCLUDED

#pragma once

#include <gsl/gsl_linalg.h>
#include<iostream>
#include<cmath>
#include<fstream>
#include <complex>
#include <ostream>
#include <math.h>
#include <vector>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_dht.h>

const double PI = 3.14159265358979323846;
const double V_Light = 299792458.;
const double eps0_vac = 8.85418781762039 * 1e-12;
const double mu0_vac = 4. * PI * 1e-7;
const std::complex<double> i(0, 1);

const std::complex<double> r1(1., 0.); ///для вакуума очев
const std::complex<double> r2(1.4504, 0.); ///SiO2 при 1мкм
const std::complex<double> r3(1.436, 9.495); ///AL при 1 мкм

const std::complex<double> eps1 = r1 * r1 * eps0_vac;
const std::complex<double> eps2 = r2 * r2 * eps0_vac;
const std::complex<double> eps3 = r3 * r3 * eps0_vac;

const double vawe_len = 1e-6; //длинна волны в вакууме
const double omega = 2. * PI / vawe_len / sqrt(eps0_vac * mu0_vac);

const std::complex<double> k_1 = omega * sqrt(eps0_vac * mu0_vac);
const std::complex<double> k_2 = k_1 * r2;
const std::complex<double> k_3 = k_1 * r3;

const std::complex<double> k_12 = k_1 * k_1;
const std::complex<double> k_22 = k_2 * k_2;
const std::complex<double> k_32 = k_3 * k_3;

const std::complex<double> sigG = std::complex<double>(1.0E-8, 1.0E-5) / eps0_vac;


const double EPS = 0.00001;




class Count_Model {
public:
	const double z_0; //Позиция источник
	const double h; //толщина слоя подложки
    std::complex<double> Fug_val_in_k1, Fug_val_in_k2, Fwg_val_in_k1, Fwg_val_in_k2;
    double steps_near_spots = 0;
    int  number_of_steps_near_spots = 0;
    std::vector<int> numbers_of_steps_in_zones_of_hankel_transforms;
    std::vector<double> steps_in_zones_of_hankel_transforms;
	int n_of_spots = 0;
	double lamda_max = 0;
	double* lamda_spots_j0 = nullptr;
	double* lamda_spots_j1 = nullptr;
	gsl_dht* dht_j0;
	gsl_dht* dht_j1;

	Count_Model(double pos, double H = 0.1);
	Count_Model(double pos, double H,
                int n_spots_near, double steps_near,
                std::vector<int> ns_spots_zones, std::vector<double> steps_in_zones);
	Count_Model(double pos, double H, 
				int n_of_spots, double lambda_max);

	std::complex<double> Count_Hankel_Tranform(std::complex<double>(*func)
                                                 (const Count_Model&, const double &, const double &, const double &), ///преобразуемая функция,
		const double & z, const double &R, double Bessel_num) const;
    std::complex<double> Count_Hankel_Tranform_with_spetial_points(std::complex<double>(*func) (const Count_Model&, const double &,
                                            const double &, const double &, const double &),
                                            const double &, const double &) const;

	std::complex<double> Count_Hankel_Tranform_with_gsl_spline(std::complex<double>(*func) (const Count_Model&, const double&,
		const double&, const double&, const double&),
		const double&, const double&, int type) const;


    std::complex<double> Count_Hankel_Tranform_near_special_spot(std::complex<double>(*func) (const Count_Model&, const double &,
                                            const double &, const double &, const double &),
                                            const double &, const double &) const;


	friend std::complex<double> U1(const Count_Model&, const double & ,const double &,const double &);
	friend std::complex<double> U0(const Count_Model&, double, const double, const double);
	friend std::complex<double> U(const Count_Model&, const double &, const double &, const double &);
	friend std::complex<double> W(const Count_Model&, const double &, const double &, const double &);
	friend std::complex<double> WP(const Count_Model&, const double &, const double &, const double &);
	friend std::complex<double> WP1(const Count_Model&, const double &, const double &, const double &);
	friend std::complex<double> WP_1(const Count_Model&, const double &, const double &, const double &);

	friend std::complex<double> d2zrF(const Count_Model&, const double & rho, const double & z, const double & z0);
	friend std::complex<double> d2zzF(const Count_Model&, const double & rho, const double & z, const double & z0);
	friend std::complex<double> dzF(const Count_Model&, const double & rho, const double & z, const double & z0);
	friend std::complex<double> drF(const Count_Model&, const double & rho, const double & z, const double & z0);
	friend std::complex<double> F(const Count_Model&, const double & rho, const double & z, const double & z0);

	friend std::complex<double> E_rho(const Count_Model&, const double &, const double &, const double &);
	friend std::complex<double> E_phi(const Count_Model&, const double &, const double &, const double &);
	friend std::complex<double> E_z(const Count_Model&, const double &, const double &, const double &);

	friend std::complex<double> H_rho(const Count_Model&, const double &, const double &, const double &);
	friend std::complex<double> H_phi(const Count_Model&, const double &, const double &, const double &);
	friend std::complex<double> H_z(const Count_Model&, const double &, const double &, const double &);

	friend std::complex<double> F_wg(const Count_Model&, const double &, const double &);
    friend std::complex<double> F_ug(const Count_Model&, const double &, const double &);

    friend std::complex<double> F_ug_n_2_is_zero(const Count_Model&, const double &, const double &);
    friend std::complex<double> F_ug_n_1_is_zero(const Count_Model&, const double &, const double &);

    friend std::complex<double> F_wg_n_2_is_zero(const Count_Model&, const double &, const double &);
    friend std::complex<double> F_wg_n_1_is_zero(const Count_Model&, const double &, const double &);

	friend std::ostream & operator << (std::ostream & out, const Count_Model&);



};
