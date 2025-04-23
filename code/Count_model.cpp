#include "Count_model.h"
#include<iomanip>



Count_Model::Count_Model(double pos, double H) : z_0(vawe_len* pos), h(vawe_len * H) {
    Fug_val_in_k1 = F_ug_n_1_is_zero(*this, k_1.real(), z_0);
    Fug_val_in_k2 = F_ug_n_2_is_zero(*this, k_2.real(), z_0);

    Fwg_val_in_k1 = F_wg_n_1_is_zero(*this, k_1.real(), z_0);
    Fwg_val_in_k2 = F_wg_n_2_is_zero(*this, k_2.real(), z_0);

    number_of_steps_near_spots = 1000;
    steps_near_spots = ((double)EPS) / number_of_steps_near_spots;

    numbers_of_steps_in_zones_of_hankel_transforms = {100000, 10000, 10000, 10000, 1000, 1000};
    steps_in_zones_of_hankel_transforms = {0.01, 1, 10, 100, 1000, 10000};

};


Count_Model::Count_Model(double pos, double H,
                         int n_spots_near, double steps_near,
                         std::vector<int> ns_spots_zones, std::vector<double> steps_in_zones) : z_0(vawe_len* pos), h(vawe_len * H) {
    Fug_val_in_k1 = F_ug_n_1_is_zero(*this, k_1.real(), z_0);
    Fug_val_in_k2 = F_ug_n_2_is_zero(*this, k_2.real(), z_0);

    Fwg_val_in_k1 = F_wg_n_1_is_zero(*this, k_1.real(), z_0);
    Fwg_val_in_k2 = F_wg_n_2_is_zero(*this, k_2.real(), z_0);

    number_of_steps_near_spots = {1000};
    steps_near_spots = {0.001};

    numbers_of_steps_in_zones_of_hankel_transforms = {100000, 10000, 10000, 10000, 1000, 1000, 100};
    steps_in_zones_of_hankel_transforms = {0.01, 1, 10, 100, 1000, 10000, 1e6};

};


Count_Model::Count_Model(double pos, double H,
	int n_spots, double lambda_max) : z_0(vawe_len* pos), h(vawe_len* H) {
	n_of_spots = n_spots;
	this->lamda_max = lamda_max;
	lamda_spots_j0 = new double[n_of_spots];
	lamda_spots_j1 = new double[n_of_spots];
	dht_j1 = gsl_dht_new(n_of_spots, 1, lambda_max);
	dht_j0 = gsl_dht_new(n_of_spots, 0, lambda_max);
	for (int i = 0; i < n_spots; ++i) {
		lamda_spots_j0[i] = gsl_dht_x_sample(dht_j0, i);
		lamda_spots_j1[i] = gsl_dht_x_sample(dht_j1, i);
	}

	Fug_val_in_k1 = F_ug_n_1_is_zero(*this, k_1.real(), z_0);
	Fug_val_in_k2 = F_ug_n_2_is_zero(*this, k_2.real(), z_0);

	Fwg_val_in_k1 = F_wg_n_1_is_zero(*this, k_1.real(), z_0);
	Fwg_val_in_k2 = F_wg_n_2_is_zero(*this, k_2.real(), z_0);

}





std::complex<double> U(const Count_Model& model,const double & z,const double & z0,const double & lambda) {
	std::complex<double> n_1 = sqrt(lambda * lambda - k_1 * k_1);
	std::complex<double> n_2 = sqrt(lambda * lambda - k_2 * k_2);
	std::complex<double> n_3 = sqrt(lambda * lambda - k_3 * k_3);

	std::complex<double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(double)2 * n_2 * model.h);
	std::complex<double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(double)2 * n_2 * model.h);

	std::complex<double> Fug = (n_2 * (R - (double)1) * (r + (double)1))
		/
		((eps2 * n_1 * (R + (double)1) - eps1 * n_2 * (R - (double)1)) * (n_2 * ((double)1 - r) + n_1 * ((double)1 + r)) - mu0_vac * sigG * sigG * n_1 * n_2 * (R - (double)1) * (r + (double)1))
		*
		(double)2 * std::complex<double>(0, 1)
		/
		(omega)
		*
		sigG;

	return Fug * exp((double) (-1) * n_1 * (z + z0));
}

std::complex<double> U1(const Count_Model& model,const double & z,const double & z0,const double & lambda) {
	std::complex<double> n_1 = sqrt(lambda * lambda - k_1 * k_1);
	std::complex<double> n_2 = sqrt(lambda * lambda - k_2 * k_2);
	std::complex<double> n_3 = sqrt(lambda * lambda - k_3 * k_3);

	std::complex<double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(double)2 * n_2 * model.h);
	std::complex<double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(double)2 * n_2 * model.h);


	return 1;
	///return Fug * exp(-n_1 * (z + z0)) * n_1;
}


std::complex<double> U0(const Count_Model& model,const double & z,const double & z0,const double & lambda) {
	std::complex<double> n_1 = sqrt(lambda * lambda - k_1 * k_1);
	std::complex<double> n_2 = sqrt(lambda * lambda - k_2 * k_2);
	std::complex<double> n_3 = sqrt(lambda * lambda - k_3 * k_3);

	std::complex<double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(double)2 * n_2 * model.h);
	std::complex<double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(double)2 * n_2 * model.h);

	std::complex<double> Fug = (n_2 * (R - (double)1) * (r + (double)1))
		/
		((eps2 * n_1 * (R + (double)1) - eps1 * n_2 * (R - (double)1)) * (n_2 * ((double)1 - r) + n_1 * ((double)1 + r)) -
   mu0_vac * sigG * sigG * n_1 * n_2 * (R - (double)1) * (r + (double)1))
		*
		(double)2 * std::complex<double>(0, 1)
		/
		(omega)
		*
		sigG;

	return Fug * exp(-n_1 * (z + z0));
}



std::complex<double> W(const Count_Model& model,const double & z,const double & z0,const double & lambda) {
	std::complex<double> n_1 = sqrt(lambda * lambda - k_1 * k_1);
	std::complex<double> n_2 = sqrt(lambda * lambda - k_2 * k_2);
	std::complex<double> n_3 = sqrt(lambda * lambda - k_3 * k_3);

	std::complex<double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(double)2 * n_2 * model.h);
	std::complex<double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(double)2 * n_2 * model.h);

	std::complex<double> Rp = R + (double)1;
	std::complex<double> Rn = R - (double)1;
	std::complex<double> rp = r + (double)1;
	std::complex<double> rn = r - (double)1;


	std::complex<double> Fwg = (double)2
		*
		(eps2 * Rp * (n_2 * (-rn) - n_1 * rp) - mu0_vac * sigG * sigG * n_2 * Rn * rp)
		/
		(((n_1 * Rp - eps1 * n_2 * Rn) * (n_2 * (-rn) + n_1 * (rp))) - mu0_vac * sigG * sigG * n_1 * n_2 * Rn * rp);

	return exp(-n_1 * abs(z - z0)) / n_1 - exp(-n_1 * abs(z + z0)) / n_1 + Fwg * exp(-n_1 * abs(z + z0)) / n_1;

}



std::complex<double> WP(const Count_Model& model,const double & z,const double & z0,const double & lambda) {
	std::complex<double> n_1 = sqrt(lambda * lambda - k_1 * k_1);
	std::complex<double> n_2 = sqrt(lambda * lambda - k_2 * k_2);
	std::complex<double> n_3 = sqrt(lambda * lambda - k_3 * k_3);

	std::complex<double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(double)2 * n_2 * model.h);
	std::complex<double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(double)2 * n_2 * model.h);

	std::complex<double> Rp = R + (double)1;
	std::complex<double> Rn = R - (double)1;
	std::complex<double> rp = r + (double)1;
	std::complex<double> rn = r - (double)1;


	std::complex<double> Fwg = (double)2
		*
		(eps2 * Rp * (n_2 * (-rn) - n_1 * rp) - mu0_vac * sigG * sigG * n_2 * Rn * rp)
		/
		(((n_1 * Rp - eps1 * n_2 * Rn) * (n_2 * (-rn) + n_1 * (rp))) - mu0_vac * sigG * sigG * n_1 * n_2 * Rn * rp);

	return Fwg * exp(-n_1 * abs(z + z0));

}


std::complex<double> WP1(const Count_Model& model,const double & z,const double & z0,const double & lambda) {
	std::complex<double> n_1 = sqrt(lambda * lambda - k_1 * k_1);
	std::complex<double> n_2 = sqrt(lambda * lambda - k_2 * k_2);
	std::complex<double> n_3 = sqrt(lambda * lambda - k_3 * k_3);

	std::complex<double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(double)2 * n_2 * model.h);
	std::complex<double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(double)2 * n_2 * model.h);

	std::complex<double> Rp = R + (double)1;
	std::complex<double> Rn = R - (double)1;
	std::complex<double> rp = r + (double)1;
	std::complex<double> rn = r - (double)1;


	std::complex<double> Fwg = (double)2
		*
		(eps2 * Rp * (n_2 * (-rn) - n_1 * rp) - mu0_vac * sigG * sigG * n_2 * Rn * rp)
		/
		(((n_1 * Rp - eps1 * n_2 * Rn) * (n_2 * (-rn) + n_1 * (rp))) - mu0_vac * sigG * sigG * n_1 * n_2 * Rn * rp);

	return Fwg * exp(-n_1 * abs(z + z0)) * n_1;

}


std::complex<double> WP_1(const Count_Model& model,const double & z,const double & z0,const double & lambda) {
	std::complex<double> n_1 = sqrt(lambda * lambda - k_1 * k_1);
	std::complex<double> n_2 = sqrt(lambda * lambda - k_2 * k_2);
	std::complex<double> n_3 = sqrt(lambda * lambda - k_3 * k_3);

	std::complex<double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(double)2 * n_2 * model.h);
	std::complex<double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(double)2 * n_2 * model.h);

	std::complex<double> Rp = R + (double)1;
	std::complex<double> Rn = R - (double)1;
	std::complex<double> rp = r + (double)1;
	std::complex<double> rn = r - (double)1;


	std::complex<double> Fwg = (double)2
		*
		(eps2 * Rp * (n_2 * (-rn) - n_1 * rp) - mu0_vac * sigG * sigG * n_2 * Rn * rp)
		/
		(((n_1 * Rp - eps1 * n_2 * Rn) * (n_2 * (-rn) + n_1 * (rp))) - mu0_vac * sigG * sigG * n_1 * n_2 * Rn * rp);

	return Fwg * exp(-n_1 * abs(z + z0)) / n_1;

}



std::complex<double> Count_Model::Count_Hankel_Tranform(std::complex<double>(*func)(const Count_Model&,const double &,const double &, const double &),
	const double & z,const double & R, double Bessel_num = 0) const {
	double lambda = 0.1;
	std::complex<double> res = 0;
	double rho = sqrt((z - z_0) * (z - z_0) + R * R);
	int k = 0;
	double step = 0.1;
	for (; lambda <= 1e9; lambda += step, ++k) {
		if (abs(lambda - k_1) < 10 || abs(lambda - k_2) < 10) continue;
		double xp = lambda - step;
		double xc = lambda;
		std::complex<double> yp = func(*this, z, z_0, xp) * (double) gsl_sf_bessel_Jnu(Bessel_num, xp * rho) * xp;
		if (Bessel_num == 1.0) {
			yp *= xp;
		}
		std::complex<double> yc = func(*this, z, z_0, xc) * (double)gsl_sf_bessel_Jnu(Bessel_num, xc * rho) * xc;
		if (Bessel_num == 1.0) {
			yc *= xc;
		}

		res += (yp + yc) / (double)2.0 * step;

		if (lambda > 1e9) {
			step = 1e4;
		}
		else if (lambda > 1e8) {
			step = 1e3;
		}
		else if (lambda > 1e7) {
			step = 1e2;
		}
		else if (lambda > 1e6) {
			step = 1e1;
		}

		if (k % 100000 == 0) {
			///std::cerr << k << " " << lambda << " " << res << std::endl;
		}

	}

	return res;


}






std::complex<double> d2zrF(const Count_Model& model, const double & rho,const double & z,const double & z0) {
	double r = sqrt(rho * rho + (z + z0) * (z + z0));
	std::complex<double> _exp = std::exp(-std::complex<double>(0, 1) * k_1 * r);
	return -rho * (z + z0) * _exp
		*
		(k_12 * (r * r) - (double)3 * std::complex<double>(0, 1) * k_1 * r - (double)3)
		/
		(r * r * r * r * r);

}



std::complex<double> d2zzF(const Count_Model& model, const double & rho,const double & z,const double & z0) {
	double r = sqrt(rho * rho + (z + z0) * (z + z0));
	std::complex<double> _exp = std::exp(-std::complex<double>(0, 1)* k_1 * r);
	std::complex<double> i(0, 1);
	return _exp / (r * r * r * r * r)
		*
		(r * r * (-k_12 * (z + z0) * (z + z0) - i * k_1 * r - (double)1) - (z + z0) * (z + z0) * (k_12 * (z + z0) * (z + z0) -
			(double)2 * i * k_1 * r - (double)2));
}


std::complex<double> dzF(const Count_Model&, const double & rho, const double & z, const double & z0) {
	double r = sqrt(rho * rho + (z + z0) * (z + z0));
	std::complex<double> _exp = std::exp(-std::complex<double>(0, 1) * k_1 * r);
	std::complex<double> i(0, 1);
	return -i * (z - z0) * _exp
		*
		(k_1 * r - i)
		/
		(r * r * r);
}


std::complex<double> drF(const Count_Model&, const double & rho, const double & z, const double & z0) {
	double r = sqrt(rho * rho + (z + z0) * (z + z0));
	std::complex<double> _exp = std::exp(- std::complex<double>(0, 1) * k_1 * r);
	std::complex<double> i(0, 1);
	return -i * rho * _exp
		*
		(k_1 * r - i)
		/
		(r * r * r);
}


std::complex<double> F(const Count_Model&, const double & rho, const double & z, const double & z0) {
	double r = sqrt(rho * rho + (z + z0) * (z + z0));
	return exp(-std::complex<double>(0, 1) * k_1 * r)
		/
		r;
}



std::complex<double> E_rho(const Count_Model& model, const double & rho, const double & z, const double & z_0) {
	std::complex<double> Koeff = (double)1
		/
		(std::complex<double>(0, 1) * omega * eps1 * mu0_vac);
	return -Koeff * d2zrF(model, rho, z, z_0) - Koeff * model.Count_Hankel_Tranform(WP, z, rho, 1.);
}





std::complex<double> E_z(const Count_Model& model, const double & rho, const double & z, const double & z_0) {
	std::complex<double> Koeff = (double)1
		/
		(std::complex<double>(0, 1) * omega * eps1 * mu0_vac);
	return -std::complex<double>(0, 1) * omega * F(model, rho, z, z_0)
		-
		Koeff * d2zzF(model, rho, z, z_0)
		+
		i * omega * model.Count_Hankel_Tranform(WP_1, z, rho, 0.)
		-
		Koeff * model.Count_Hankel_Tranform(WP1, z, rho, 0.);

}


std::complex<double> H_rho(const Count_Model& model, const double & rho, const double & z, const double & z_0) {
	return -model.Count_Hankel_Tranform(U0, z, rho, 1.)
		/
		mu0_vac;
}



std::complex<double>H_phi (const Count_Model& model, const double & rho, const double & z, const double & z_0) {
	return -drF(model, rho, z, z_0) / mu0_vac
		+
		model.Count_Hankel_Tranform(WP_1, z, rho, 1.) / mu0_vac;
}


std::complex<double>H_z(const Count_Model& model, const double & rho, const double & z, const double & z_0) {
	return k_12 / mu0_vac * model.Count_Hankel_Tranform(U, z, rho, 0.)
		+
		k_12 / mu0_vac + model.Count_Hankel_Tranform(U1, z, rho, 0);
}



std::complex<double> F_wg(const Count_Model& model, const double & lambda, const double & z0) {
    std::complex<double> n_1 = sqrt(lambda * lambda - k_1 * k_1);
	std::complex<double> n_2 = sqrt(lambda * lambda - k_2 * k_2);
	std::complex<double> n_3 = sqrt(lambda * lambda - k_3 * k_3);

	std::complex<double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(double)2 * n_2 * model.h);
	std::complex<double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(double)2 * n_2 * model.h);

	std::complex<double> Rp = R + (double)1;
	std::complex<double> Rn = R - (double)1;
	std::complex<double> rp = r + (double)1;
	std::complex<double> rn = r - (double)1;


	std::complex<double> Fwg = (double)2
		*
		(eps2 * Rp * (n_2 * (-rn) - n_1 * rp) - mu0_vac * sigG * sigG * n_2 * Rn * rp)
		/
		(((n_1 * Rp - eps1 * n_2 * Rn) * (n_2 * (-rn) + n_1 * (rp))) - mu0_vac * sigG * sigG * n_1 * n_2 * Rn * rp);

    return Fwg;
}



std::complex<double> F_ug(const Count_Model& model, const double & lambda, const double & z0) {
    ///std::cout << "lambda = " << lambda << '\n';
    std::complex<double> n_1 = sqrt(lambda * lambda - k_1 * k_1);
	std::complex<double> n_2 = sqrt(lambda * lambda - k_2 * k_2);
	std::complex<double> n_3 = sqrt(lambda * lambda - k_3 * k_3);

	std::complex<double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(double)2 * n_2 * model.h);
	std::complex<double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(double)2 * n_2 * model.h);

	//std::cout << std::fixed << std::setprecision(40);
	//std::cout << "lambda: " << lambda << '\n';
	//std::cout << "k_1:" << k_1 << '\n';
	//std::cout << "lambda^2: " << lambda * lambda << '\n';
	//std::cout << "k_1^2:" << k_12 << '\n';
	//std::cout << "n_1: " << n_1 << '\n';
	//std::cout << "n_1^2: " << n_1*n_1 << '\n';




    std::complex<double> Fug = (n_2 * (R - (double)1) * (r + (double)1))
		/
		( (eps2 * n_1 * (R + (double)1) - eps1 * n_2 * (R - (double)1)) * (n_2 * ((double)1 - r) + n_1 * ((double)1 + r)) -
        (mu0_vac * sigG * sigG * n_1 * n_2 * (R - (double)1) * (r + (double)1)))
		*
		(double)2 * i
		/
		(omega)
		*
		sigG;

    if (abs(n_1) <= EPS) {
        //std::cout << "its_me " << lambda << ' ' << model.steps_near_spots << '\n';
        return F_ug_n_1_is_zero(model, lambda, z0);
    }
    if (abs(n_2) <= EPS) {
        ///костыль
        //std::cout << "its_me " << lambda << '\n';
        ///return {0, 0};
        return F_ug_n_2_is_zero(model, lambda, z0);
    }
    return Fug;
}

std::complex<double> F_ug_n_2_is_zero(const Count_Model& model, const double & lambda, const double & z0) {
    std::complex<double> n_1 = sqrt(lambda * lambda - k_1 * k_1);
	std::complex<double> n_2 = sqrt(lambda * lambda - k_2 * k_2);
	std::complex<double> n_3 = sqrt(lambda * lambda - k_3 * k_3);
    std::complex<double> C_eps = model.h + eps3 / eps2 / n_3;
    std::complex<double> C = model.h + (double) 1 / n_3;

    return (double) (-4) * C * i * sigG / omega /
    (((double) 2 * eps2 * n_1 * C_eps + eps1) * ((double) 2) *((double) 1 + n_2) +
     (double) 4 * mu0_vac * sigG * sigG * n_1 * C);

}


std::complex<double> F_ug_n_1_is_zero(const Count_Model& model, const double & lambda, const double & z0) {
    std::complex<double> n_1 = sqrt(lambda * lambda - k_1 * k_1);
	std::complex<double> n_2 = sqrt(lambda * lambda - k_2 * k_2);
	std::complex<double> n_3 = sqrt(lambda * lambda - k_3 * k_3);
    std::complex<double> C_eps = model.h + eps3 / eps2 / n_3;
    std::complex<double> C = model.h + (double) 1 / n_3;

    std::complex<double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(double)2 * n_2 * model.h);
	std::complex<double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(double)2 * n_2 * model.h);

    return (double) 2 * i * sigG / omega *
    (r + (double) 1) /
    (-eps1 * n_2 * ((double) 1 - r))  ;

}


std::complex<double> F_wg_n_1_is_zero(const Count_Model& model, const double & lambda, const double & z0) {
    std::complex<double> n_1 = sqrt(lambda * lambda - k_1 * k_1);
	std::complex<double> n_2 = sqrt(lambda * lambda - k_2 * k_2);
	std::complex<double> n_3 = sqrt(lambda * lambda - k_3 * k_3);
    std::complex<double> C_eps = model.h + eps3 / eps2 / n_3;
    std::complex<double> C = model.h + (double) 1 / n_3;

    std::complex<double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(double)2 * n_2 * model.h);
	std::complex<double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(double)2 * n_2 * model.h);

    return  (double) 2 * (mu0_vac * sigG * sigG * (r + (double) 1)
            /
            (eps1 * n_2 * ((double) 1 - r))
            +
            eps2 * (r + (double) 1)
            /
            (eps1 * (R - (double) 1) * n_2));
}





std::complex<double> F_wg_n_2_is_zero(const Count_Model& model, const double & lambda, const double & z0) {
    std::complex<double> n_1 = sqrt(lambda * lambda - k_1 * k_1);
	std::complex<double> n_2 = sqrt(lambda * lambda - k_2 * k_2);
	std::complex<double> n_3 = sqrt(lambda * lambda - k_3 * k_3);
    std::complex<double> C_eps = model.h + eps3 / eps2 / n_3;
    std::complex<double> C = model.h + (double) 1 / n_3;

    std::complex<double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(double)2 * n_2 * model.h);
	std::complex<double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(double)2 * n_2 * model.h);

    return (double) 2 * (eps2 * C_eps * ((double) 1 + n_1 * C) + mu0_vac * sigG * sigG * C)
            /
            (mu0_vac * sigG * sigG * n_1 * C - (eps2 * n_1 * C_eps + eps1) * ((double) 1 + n_1 * C));
}


std::complex<double> Underintegral_d2Az_dzdr(const Count_Model& model, const double & lambda,
                                                  const double & z, const double & z0, const double & rho) {
    std::complex<double> n_1 = sqrt(lambda * lambda - k_1 * k_1);
    double R = sqrt((z - z0) * (z - z0) + rho * rho);

    ///std::cerr << "z: " << z << '\n';
    ///std::cerr << "n_1: " << n_1 << '\n';
	std::cerr << R;


	if (lambda * rho == 0.0) {
		return (double)gsl_sf_bessel_Jnu(1., 0.000000001)* F_ug(model, lambda, z0)* exp(-n_1 * (z + z0))* lambda* lambda;
	}
    return  (double)gsl_sf_bessel_Jnu(1., lambda * rho) * F_ug(model, lambda, z0) * exp(-n_1 * (z + z0)) * lambda * lambda;
}



std::complex<double> E_phi(const Count_Model& model,
                                const double & z, const double & z0, const double & rho) {
    return model.Count_Hankel_Tranform_with_gsl_spline(Underintegral_d2Az_dzdr, z, rho, 1);
}


//Преобразования Ханкеля с учетом двух "особенных отчек"
std::complex<double> Count_Model::Count_Hankel_Tranform_with_spetial_points
                                         (std::complex<double>(*func) (const Count_Model&, const double &,
                                            const double &, const double &, const double &),
            const double & z, const double & R) const {

	std::cout << z << std::endl;
    //отдельно обсчитываем особенность
    std::complex<double> val_near_k1 = {0, 0};
    std::complex<double> val_near_k2 = {0, 0};

    for (int i = -number_of_steps_near_spots; i < number_of_steps_near_spots; ++i) {
        if (i == 0) continue;
        //std::cout << steps_near_spots << '\n';
        val_near_k1 += func(*this, k_1.real() + i * steps_near_spots, z, this->z_0, R) * steps_near_spots;
        val_near_k2 += func(*this, k_2.real() + i * steps_near_spots, z, this->z_0, R) * steps_near_spots;
    }

    val_near_k1 += func(*this, k_1.real(), z, this->z_0, R) * steps_near_spots;
    val_near_k2 += func(*this, k_2.real(), z, this->z_0, R) * steps_near_spots;
    //std::cout << val_near_k1 << ' ' << val_near_k2 << std::endl;


    std::complex<double> res = {0, 0};

    double lambda = 0.0000000001;
    for (int zone_num = 0; zone_num < numbers_of_steps_in_zones_of_hankel_transforms.size(); ++zone_num) {
        for (int i = 0; i < numbers_of_steps_in_zones_of_hankel_transforms[zone_num]; ++i) {
            if (lambda >= k_1.real() - number_of_steps_near_spots * steps_near_spots &&
                lambda <= k_1.real() + number_of_steps_near_spots * steps_near_spots)
            {
                lambda += steps_in_zones_of_hankel_transforms[zone_num];
                continue;
            }
            if (lambda >= k_2.real() - number_of_steps_near_spots * steps_near_spots &&
                lambda <= k_2.real() + number_of_steps_near_spots * steps_near_spots)
            {
                lambda += steps_in_zones_of_hankel_transforms[zone_num];
                continue;
            }
            res += func(*this, lambda, z, this->z_0, R) * steps_in_zones_of_hankel_transforms[zone_num];
            lambda += steps_in_zones_of_hankel_transforms[zone_num];
        }
    }

    //std::cout << '\n' << res << '\n';

    return res + val_near_k1 + val_near_k2;
}


std::complex<double> Count_Model::Count_Hankel_Tranform_with_gsl_spline
(std::complex<double>(*func) (const Count_Model&, const double&,
	const double&, const double&, const double&),
	const double& z, const double& R, int type) const {
	std::complex<double> res1 = 0;
	std::complex<double> res0 = 0;
	for (int i = 0; i < n_of_spots - 1; ++i) {
		res1 += func(*this, lamda_spots_j1[i], z, this->z_0, R) * (lamda_spots_j1[i + 1] - lamda_spots_j1[i]);
		res0 += func(*this, lamda_spots_j0[i], z, this->z_0, R) * (lamda_spots_j0[i + 1] - lamda_spots_j0[i]);
	}

	if (type) {
		return res1;
	}
	else {
		return res0;
	}
}



std::ostream & operator << (std::ostream & out, const std::vector<double> & v) {
    for (auto el : v) {
        out << el << ' ';
    }

    return out;
}



std::ostream & operator << (std::ostream & out, const Count_Model& model) {
    out << "\n|----------------------------- \n| Model info \n";
    out << "| params of layers: \n";
    out << "| z0: " << model.z_0 << "\n";
    out << "| H: " << model.h << "\n";
    out << "| \n";
    out << "| k_1: " << k_1 << std::endl;
    out << "| Fug in k1: " << model.Fug_val_in_k1 << std::endl;
    out << "| Fwg in k1: " << model.Fwg_val_in_k1 << std::endl;
    out << "| k_2: " << k_2 << std::endl;
    out << "| Fug in k2: " << model.Fug_val_in_k2 << std::endl;
    out << "| Fwg in k2: " << model.Fwg_val_in_k2 << std::endl;
    out << "| k_3: " << k_3 << std::endl;
    out << "|\n| Model params\n";
    out << "| step value near spots: " << model.steps_near_spots << "\n";
    out << "| number of steps near spots: " << model.number_of_steps_near_spots << "\n";
    out << "| number of counting zones: " << model.steps_in_zones_of_hankel_transforms.size() << "\n";
    out << "| zones info:\n";
    for (int i = 0; i < model.steps_in_zones_of_hankel_transforms.size(); ++i) {
        out << "|| zone" << i <<  ": n spots = " << model.numbers_of_steps_in_zones_of_hankel_transforms[i] << ' ';
        out << " step size = " << model.steps_in_zones_of_hankel_transforms[i] << "\n";
    }
    out << "|----------------------------- \n \n";

    return out;
}




