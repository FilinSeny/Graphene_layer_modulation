#include "Count_model.h"

Count_Model::Count_Model(long double pos, long double H) : z_0(vawe_len* pos), h(vawe_len * H) {
    Fug_val_in_k1 = F_ug_n_1_is_zero(*this, k_1.real(), z_0);
    Fug_val_in_k2 = F_ug_n_2_is_zero(*this, k_2.real(), z_0);

    Fwg_val_in_k1 = F_wg_n_1_is_zero(*this, k_1.real(), z_0);
    Fwg_val_in_k2 = F_wg_n_2_is_zero(*this, k_2.real(), z_0);

    number_of_steps_near_spots = 1000;
    steps_near_spots = 0.001;

    numbers_of_steps_in_zones_of_hankel_transforms = {100000, 10000, 10000, 10000, 1000, 1000};
    steps_in_zones_of_hankel_transforms = {0.01, 1, 10, 100, 1000, 10000};

};


Count_Model::Count_Model(long double pos, long double H,
                         int n_spots_near, long double steps_near,
                         std::vector<int> ns_spots_zones, std::vector<long double> steps_in_zones) : z_0(vawe_len* pos), h(vawe_len * H) {
    Fug_val_in_k1 = F_ug_n_1_is_zero(*this, k_1.real(), z_0);
    Fug_val_in_k2 = F_ug_n_2_is_zero(*this, k_2.real(), z_0);

    Fwg_val_in_k1 = F_wg_n_1_is_zero(*this, k_1.real(), z_0);
    Fwg_val_in_k2 = F_wg_n_2_is_zero(*this, k_2.real(), z_0);

    number_of_steps_near_spots = {1000};
    steps_near_spots = {0.001};

    numbers_of_steps_in_zones_of_hankel_transforms = {100000, 10000, 10000, 10000, 1000, 1000, 100};
    steps_in_zones_of_hankel_transforms = {0.01, 1, 10, 100, 1000, 10000, 1e6};

};





std::complex<long double> U(const Count_Model& model,const long double & z,const long double & z0,const long double & lambda) {
	std::complex<long double> n_1 = sqrt(lambda * lambda - k_1 * k_1);
	std::complex<long double> n_2 = sqrt(lambda * lambda - k_2 * k_2);
	std::complex<long double> n_3 = sqrt(lambda * lambda - k_3 * k_3);

	std::complex<long double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(long double)2 * n_2 * model.h);
	std::complex<long double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(long double)2 * n_2 * model.h);

	std::complex<long double> Fug = (n_2 * (R - (long double)1) * (r + (long double)1))
		/
		((eps2 * n_1 * (R + (long double)1) - eps1 * n_2 * (R - (long double)1)) * (n_2 * ((long double)1 - r) + n_1 * ((long double)1 + r)) - mu0_vac * sigG * sigG * n_1 * n_2 * (R - (long double)1) * (r + (long double)1))
		*
		(long double)2 * std::complex<long double>(0, 1)
		/
		(omega)
		*
		sigG;

	return Fug * exp((long double) (-1) * n_1 * (z + z0));
}

std::complex<long double> U1(const Count_Model& model,const long double & z,const long double & z0,const long double & lambda) {
	std::complex<long double> n_1 = sqrt(lambda * lambda - k_1 * k_1);
	std::complex<long double> n_2 = sqrt(lambda * lambda - k_2 * k_2);
	std::complex<long double> n_3 = sqrt(lambda * lambda - k_3 * k_3);

	std::complex<long double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(long double)2 * n_2 * model.h);
	std::complex<long double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(long double)2 * n_2 * model.h);

	std::complex<long double> Fug = (n_2 * (R - (long double)1) * (r + (long double)1))
		/
		((eps2 * n_1 * (R + (long double)1) - eps1 * n_2 * (R - (long double)1)) * (n_2 * ((long double)1 - r) + n_1 * ((long double)1 + r)) - mu0_vac * sigG * sigG * n_1 * n_2 * (R - (long double)1) * (r + (long double)1))
		*
		///��� ������!
		(long double)2 * std::complex<long double>(0, 1)
		/
		(omega)
		*
		sigG;

	return Fug * exp(-n_1 * (z + z0)) * n_1;
}


std::complex<long double> U0(const Count_Model& model,const long double & z,const long double & z0,const long double & lambda) {
	std::complex<long double> n_1 = sqrt(lambda * lambda - k_1 * k_1);
	std::complex<long double> n_2 = sqrt(lambda * lambda - k_2 * k_2);
	std::complex<long double> n_3 = sqrt(lambda * lambda - k_3 * k_3);

	std::complex<long double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(long double)2 * n_2 * model.h);
	std::complex<long double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(long double)2 * n_2 * model.h);

	std::complex<long double> Fug = (n_2 * (R - (long double)1) * (r + (long double)1))
		/
		((eps2 * n_1 * (R + (long double)1) - eps1 * n_2 * (R - (long double)1)) * (n_2 * ((long double)1 - r) + n_1 * ((long double)1 + r)) -
   mu0_vac * sigG * sigG * n_1 * n_2 * (R - (long double)1) * (r + (long double)1))
		*
		(long double)2 * std::complex<long double>(0, 1)
		/
		(omega)
		*
		sigG;

	return Fug * exp(-n_1 * (z + z0));
}



std::complex<long double> W(const Count_Model& model,const long double & z,const long double & z0,const long double & lambda) {
	std::complex<long double> n_1 = sqrt(lambda * lambda - k_1 * k_1);
	std::complex<long double> n_2 = sqrt(lambda * lambda - k_2 * k_2);
	std::complex<long double> n_3 = sqrt(lambda * lambda - k_3 * k_3);

	std::complex<long double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(long double)2 * n_2 * model.h);
	std::complex<long double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(long double)2 * n_2 * model.h);

	std::complex<long double> Rp = R + (long double)1;
	std::complex<long double> Rn = R - (long double)1;
	std::complex<long double> rp = r + (long double)1;
	std::complex<long double> rn = r - (long double)1;


	std::complex<long double> Fwg = (long double)2
		*
		(eps2 * Rp * (n_2 * (-rn) - n_1 * rp) - mu0_vac * sigG * sigG * n_2 * Rn * rp)
		/
		(((n_1 * Rp - eps1 * n_2 * Rn) * (n_2 * (-rn) + n_1 * (rp))) - mu0_vac * sigG * sigG * n_1 * n_2 * Rn * rp);

	return exp(-n_1 * abs(z - z0)) / n_1 - exp(-n_1 * abs(z + z0)) / n_1 + Fwg * exp(-n_1 * abs(z + z0)) / n_1;

}



std::complex<long double> WP(const Count_Model& model,const long double & z,const long double & z0,const long double & lambda) {
	std::complex<long double> n_1 = sqrt(lambda * lambda - k_1 * k_1);
	std::complex<long double> n_2 = sqrt(lambda * lambda - k_2 * k_2);
	std::complex<long double> n_3 = sqrt(lambda * lambda - k_3 * k_3);

	std::complex<long double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(long double)2 * n_2 * model.h);
	std::complex<long double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(long double)2 * n_2 * model.h);

	std::complex<long double> Rp = R + (long double)1;
	std::complex<long double> Rn = R - (long double)1;
	std::complex<long double> rp = r + (long double)1;
	std::complex<long double> rn = r - (long double)1;


	std::complex<long double> Fwg = (long double)2
		*
		(eps2 * Rp * (n_2 * (-rn) - n_1 * rp) - mu0_vac * sigG * sigG * n_2 * Rn * rp)
		/
		(((n_1 * Rp - eps1 * n_2 * Rn) * (n_2 * (-rn) + n_1 * (rp))) - mu0_vac * sigG * sigG * n_1 * n_2 * Rn * rp);

	return Fwg * exp(-n_1 * abs(z + z0));

}


std::complex<long double> WP1(const Count_Model& model,const long double & z,const long double & z0,const long double & lambda) {
	std::complex<long double> n_1 = sqrt(lambda * lambda - k_1 * k_1);
	std::complex<long double> n_2 = sqrt(lambda * lambda - k_2 * k_2);
	std::complex<long double> n_3 = sqrt(lambda * lambda - k_3 * k_3);

	std::complex<long double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(long double)2 * n_2 * model.h);
	std::complex<long double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(long double)2 * n_2 * model.h);

	std::complex<long double> Rp = R + (long double)1;
	std::complex<long double> Rn = R - (long double)1;
	std::complex<long double> rp = r + (long double)1;
	std::complex<long double> rn = r - (long double)1;


	std::complex<long double> Fwg = (long double)2
		*
		(eps2 * Rp * (n_2 * (-rn) - n_1 * rp) - mu0_vac * sigG * sigG * n_2 * Rn * rp)
		/
		(((n_1 * Rp - eps1 * n_2 * Rn) * (n_2 * (-rn) + n_1 * (rp))) - mu0_vac * sigG * sigG * n_1 * n_2 * Rn * rp);

	return Fwg * exp(-n_1 * abs(z + z0)) * n_1;

}


std::complex<long double> WP_1(const Count_Model& model,const long double & z,const long double & z0,const long double & lambda) {
	std::complex<long double> n_1 = sqrt(lambda * lambda - k_1 * k_1);
	std::complex<long double> n_2 = sqrt(lambda * lambda - k_2 * k_2);
	std::complex<long double> n_3 = sqrt(lambda * lambda - k_3 * k_3);

	std::complex<long double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(long double)2 * n_2 * model.h);
	std::complex<long double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(long double)2 * n_2 * model.h);

	std::complex<long double> Rp = R + (long double)1;
	std::complex<long double> Rn = R - (long double)1;
	std::complex<long double> rp = r + (long double)1;
	std::complex<long double> rn = r - (long double)1;


	std::complex<long double> Fwg = (long double)2
		*
		(eps2 * Rp * (n_2 * (-rn) - n_1 * rp) - mu0_vac * sigG * sigG * n_2 * Rn * rp)
		/
		(((n_1 * Rp - eps1 * n_2 * Rn) * (n_2 * (-rn) + n_1 * (rp))) - mu0_vac * sigG * sigG * n_1 * n_2 * Rn * rp);

	return Fwg * exp(-n_1 * abs(z + z0)) / n_1;

}



std::complex<long double> Count_Model::Count_Hankel_Tranform(std::complex<long double>(*func)(const Count_Model&,const long double &,const long double &, const long double &),
	const long double & z,const long double & R, double Bessel_num = 0) const {
	long double lambda = 0.1;
	std::complex<long double> res = 0;
	long double rho = sqrt((z - z_0) * (z - z_0) + R * R);
	int k = 0;
	long double step = 0.1;
	for (; lambda <= 1e9; lambda += step, ++k) {
		if (abs(lambda - k_1) < 10 || abs(lambda - k_2) < 10) continue;
		long double xp = lambda - step;
		long double xc = lambda;
		std::complex<long double> yp = func(*this, z, z_0, xp) * std::cyl_bessel_j(Bessel_num, xp * rho) * xp;
		if (Bessel_num == 1.0) {
			yp *= xp;
		}
		std::complex<long double> yc = func(*this, z, z_0, xc) * std::cyl_bessel_j(Bessel_num, xc * rho) * xc;
		if (Bessel_num == 1.0) {
			yc *= xc;
		}

		res += (yp + yc) / (long double)2.0 * step;

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



std::complex<long double> Count_Model::Count_Hankel_Tranform_near_special_spot(std::complex<long double>(*func)
                                                                               (const Count_Model&,const long double &,const long double &, const long double &),
	const long double & z,const long double & R, std::string func_name, int number_of_point, double Bessel_num = 0) const {

    long double spot = number_of_point == 1 ? k_1.real() : k_2.real();
    std::complex<long double> spot_val {0, 0};
    if (func_name == "Fug") {
        spot_val = number_of_point == 1 ? Fug_val_in_k1 : Fug_val_in_k2;
    } else {
        spot_val = number_of_point == 1 ? Fwg_val_in_k1 : Fwg_val_in_k2;
    }

    std::complex<long double> res {0, 0};
    long double rho = sqrt((z - z_0) * (z - z_0) + R * R);

    for (int i = -number_of_steps_near_spots; i < number_of_steps_near_spots; ++i) {
        if (!i) continue;
        long double xp = spot - (steps_near_spots * (i - 1));
		long double xc = spot - (steps_near_spots * (i));
		std::complex<long double> yp = func(*this, z, z_0, xp) * std::cyl_bessel_j(Bessel_num, xp * rho) * xp;
		if (Bessel_num == 1.0) {
			yp *= xp;
		}
		std::complex<long double> yc = func(*this, z, z_0, xc) * std::cyl_bessel_j(Bessel_num, xc * rho) * xc;
		if (Bessel_num == 1.0) {
			yc *= xc;
		}

		res += (yp + yc) / (long double)2.0 * steps_near_spots;
    }

    return res + steps_near_spots * spot_val;
}



//�������������� ������� � ������ ���� "��������� �����"
std::complex<long double> Count_Model::Count_Hankel_Tranform_with_spetial_points(std::complex<long double>(*func)(const Count_Model&,const long double &,const long double &, const long double &),
	const long double & z,const long double & R, std::string func_name, double Bessel_num = 0) const {

    //�������� ����������� �����������
    std::complex<long double> val_near_k1 = Count_Hankel_Tranform_near_special_spot(func, z, R, func_name, 1, Bessel_num);
    std::complex<long double> val_near_k2 = Count_Hankel_Tranform_near_special_spot(func, z, R, func_name, 2, Bessel_num);

    long double add = 0;
    int num_of_zone = 0;
    long double lambda = 0;
    std::complex<long double> res = {0, 0};
    long double rho = sqrt((z - z_0) * (z - z_0) + R * R);

    for (; num_of_zone < numbers_of_steps_in_zones_of_hankel_transforms.size(); ++num_of_zone) {
        for (int i = 0; i < numbers_of_steps_in_zones_of_hankel_transforms[num_of_zone]; ++i) {
            if (lambda >= k_1.real() - number_of_steps_near_spots * steps_near_spots &&
                 lambda <= k_1.real() + number_of_steps_near_spots * steps_near_spots)
            {
                 continue;
            }
            if (lambda >= k_2.real() - number_of_steps_near_spots * steps_near_spots &&
                 lambda <= k_2.real() + number_of_steps_near_spots * steps_near_spots)
            {
                 continue;
            }
            res += func(*this, z, z_0, lambda) * std::cyl_bessel_j(Bessel_num, lambda * rho) * lambda *
            steps_in_zones_of_hankel_transforms[num_of_zone];
            lambda += steps_in_zones_of_hankel_transforms[num_of_zone];
        }
    }

    return res + val_near_k1 + val_near_k2;
}



std::complex<long double> d2zrF(const Count_Model& model, const long double & rho,const long double & z,const long double & z0) {
	long double r = sqrt(rho * rho + (z + z0) * (z + z0));
	std::complex<long double> _exp = std::exp(-std::complex<long double>(0, 1) * k_1 * r);
	return -rho * (z + z0) * _exp
		*
		(k_12 * (r * r) - (long double)3 * std::complex<long double>(0, 1) * k_1 * r - (long double)3)
		/
		(r * r * r * r * r);

}



std::complex<long double> d2zzF(const Count_Model& model, const long double & rho,const long double & z,const long double & z0) {
	long double r = sqrt(rho * rho + (z + z0) * (z + z0));
	std::complex<long double> _exp = std::exp(-std::complex<long double>(0, 1)* k_1 * r);
	std::complex<long double> i(0, 1);
	return _exp / (r * r * r * r * r)
		*
		(r * r * (-k_12 * (z + z0) * (z + z0) - i * k_1 * r - (long double)1) - (z + z0) * (z + z0) * (k_12 * (z + z0) * (z + z0) -
			(long double)2 * i * k_1 * r - (long double)2));
}


std::complex<long double> dzF(const Count_Model&, const long double & rho, const long double & z, const long double & z0) {
	long double r = sqrt(rho * rho + (z + z0) * (z + z0));
	std::complex<long double> _exp = std::exp(-std::complex<long double>(0, 1) * k_1 * r);
	std::complex<long double> i(0, 1);
	return -i * (z - z0) * _exp
		*
		(k_1 * r - i)
		/
		(r * r * r);
}


std::complex<long double> drF(const Count_Model&, const long double & rho, const long double & z, const long double & z0) {
	long double r = sqrt(rho * rho + (z + z0) * (z + z0));
	std::complex<long double> _exp = std::exp(- std::complex<long double>(0, 1) * k_1 * r);
	std::complex<long double> i(0, 1);
	return -i * rho * _exp
		*
		(k_1 * r - i)
		/
		(r * r * r);
}


std::complex<long double> F(const Count_Model&, const long double & rho, const long double & z, const long double & z0) {
	long double r = sqrt(rho * rho + (z + z0) * (z + z0));
	return exp(-std::complex<long double>(0, 1) * k_1 * r)
		/
		r;
}



std::complex<long double> E_rho(const Count_Model& model, const long double & rho, const long double & z, const long double & z_0) {
	std::complex<long double> Koeff = (long double)1
		/
		(std::complex<long double>(0, 1) * omega * eps1 * mu0_vac);
	return -Koeff * d2zrF(model, rho, z, z_0) - Koeff * model.Count_Hankel_Tranform(WP, z, rho, 1.);
}



std::complex<long double> E_phi(const Count_Model& model, const long double & rho, const long double & z, const long double & z_0) {
	return -std::complex<long double>(0, 1) * omega * model.Count_Hankel_Tranform(U, z, rho, 1.);
}



std::complex<long double> E_z(const Count_Model& model, const long double & rho, const long double & z, const long double & z_0) {
	std::complex<long double> Koeff = (long double)1
		/
		(std::complex<long double>(0, 1) * omega * eps1 * mu0_vac);
	return -std::complex<long double>(0, 1) * omega * F(model, rho, z, z_0)
		-
		Koeff * d2zzF(model, rho, z, z_0)
		+
		i * omega * model.Count_Hankel_Tranform(WP_1, z, rho, 0.)
		-
		Koeff * model.Count_Hankel_Tranform(WP1, z, rho, 0.);

}


std::complex<long double> H_rho(const Count_Model& model, const long double & rho, const long double & z, const long double & z_0) {
	return -model.Count_Hankel_Tranform(U0, z, rho, 1.)
		/
		mu0_vac;
}



std::complex<long double>H_phi (const Count_Model& model, const long double & rho, const long double & z, const long double & z_0) {
	return -drF(model, rho, z, z_0) / mu0_vac
		+
		model.Count_Hankel_Tranform(WP_1, z, rho, 1.) / mu0_vac;
}


std::complex<long double>H_z(const Count_Model& model, const long double & rho, const long double & z, const long double & z_0) {
	return k_12 / mu0_vac * model.Count_Hankel_Tranform(U, z, rho, 0.)
		+
		k_12 / mu0_vac + model.Count_Hankel_Tranform(U1, z, rho, 0);
}



std::complex<long double> F_wg(const Count_Model& model, const long double & lambda, const long double & z0) {
    std::complex<long double> n_1 = sqrt(lambda * lambda - k_1 * k_1);
	std::complex<long double> n_2 = sqrt(lambda * lambda - k_2 * k_2);
	std::complex<long double> n_3 = sqrt(lambda * lambda - k_3 * k_3);

	std::complex<long double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(long double)2 * n_2 * model.h);
	std::complex<long double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(long double)2 * n_2 * model.h);

	std::complex<long double> Rp = R + (long double)1;
	std::complex<long double> Rn = R - (long double)1;
	std::complex<long double> rp = r + (long double)1;
	std::complex<long double> rn = r - (long double)1;


	std::complex<long double> Fwg = (long double)2
		*
		(eps2 * Rp * (n_2 * (-rn) - n_1 * rp) - mu0_vac * sigG * sigG * n_2 * Rn * rp)
		/
		(((n_1 * Rp - eps1 * n_2 * Rn) * (n_2 * (-rn) + n_1 * (rp))) - mu0_vac * sigG * sigG * n_1 * n_2 * Rn * rp);

    return Fwg;
}



std::complex<long double> F_ug(const Count_Model& model, const long double & lambda, const long double & z0) {
    std::complex<long double> n_1 = sqrt(lambda * lambda - k_1 * k_1);
	std::complex<long double> n_2 = sqrt(lambda * lambda - k_2 * k_2);
	std::complex<long double> n_3 = sqrt(lambda * lambda - k_3 * k_3);

	std::complex<long double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(long double)2 * n_2 * model.h);
	std::complex<long double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(long double)2 * n_2 * model.h);

	///std::cout << R.real() << std::endl;
	std::complex<long double> Fug = (n_2 * (R - (long double)1) * (r + (long double)1))
		/
		( (eps2 * n_1 * (R + (long double)1) - eps1 * n_2 * (R - (long double)1)) * (n_2 * ((long double)1 - r) + n_1 * ((long double)1 + r)) -
        (mu0_vac * sigG * sigG * n_1 * n_2 * (R - (long double)1) * (r + (long double)1)))
		*
		(long double)2 * i
		/
		(omega)
		*
		sigG;
    return Fug;
}

std::complex<long double> F_ug_n_2_is_zero(const Count_Model& model, const long double & lambda, const long double & z0) {
    std::complex<long double> n_1 = sqrt(lambda * lambda - k_1 * k_1);
	std::complex<long double> n_2 = sqrt(lambda * lambda - k_2 * k_2);
	std::complex<long double> n_3 = sqrt(lambda * lambda - k_3 * k_3);
    std::complex<long double> C_eps = model.h + eps3 / eps2 / n_3;
    std::complex<long double> C = model.h + (long double) 1 / n_3;

    return (long double) (-4) * C * i * sigG / omega /
    (((long double) 2 * eps2 * n_1 * C_eps + eps1) * ((long double) 2) *((long double) 1 + n_2) +
     (long double) 4 * mu0_vac * sigG * sigG * n_1 * C);

}



std::complex<long double> F_ug_n_1_is_zero(const Count_Model& model, const long double & lambda, const long double & z0) {
    std::complex<long double> n_1 = sqrt(lambda * lambda - k_1 * k_1);
	std::complex<long double> n_2 = sqrt(lambda * lambda - k_2 * k_2);
	std::complex<long double> n_3 = sqrt(lambda * lambda - k_3 * k_3);
    std::complex<long double> C_eps = model.h + eps3 / eps2 / n_3;
    std::complex<long double> C = model.h + (long double) 1 / n_3;

    std::complex<long double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(long double)2 * n_2 * model.h);
	std::complex<long double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(long double)2 * n_2 * model.h);

    return (long double) 2 * i * sigG / omega *
    (r + (long double) 1) /
    (-eps1 * n_2 * ((long double) 1 - r))  ;

}


std::complex<long double> F_wg_n_1_is_zero(const Count_Model& model, const long double & lambda, const long double & z0) {
    std::complex<long double> n_1 = sqrt(lambda * lambda - k_1 * k_1);
	std::complex<long double> n_2 = sqrt(lambda * lambda - k_2 * k_2);
	std::complex<long double> n_3 = sqrt(lambda * lambda - k_3 * k_3);
    std::complex<long double> C_eps = model.h + eps3 / eps2 / n_3;
    std::complex<long double> C = model.h + (long double) 1 / n_3;

    std::complex<long double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(long double)2 * n_2 * model.h);
	std::complex<long double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(long double)2 * n_2 * model.h);

    return  (long double) 2 * (mu0_vac * sigG * sigG * (r + (long double) 1)
            /
            (eps1 * n_2 * ((long double) 1 - r))
            +
            eps2 * (r + (long double) 1)
            /
            (eps1 * (R - (long double) 1) * n_2));
}


std::complex<long double> F_wg_n_2_is_zero(const Count_Model& model, const long double & lambda, const long double & z0) {
    std::complex<long double> n_1 = sqrt(lambda * lambda - k_1 * k_1);
	std::complex<long double> n_2 = sqrt(lambda * lambda - k_2 * k_2);
	std::complex<long double> n_3 = sqrt(lambda * lambda - k_3 * k_3);
    std::complex<long double> C_eps = model.h + eps3 / eps2 / n_3;
    std::complex<long double> C = model.h + (long double) 1 / n_3;

    std::complex<long double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(long double)2 * n_2 * model.h);
	std::complex<long double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(long double)2 * n_2 * model.h);

    return (long double) 2 * (eps2 * C_eps * ((long double) 1 + n_1 * C) + mu0_vac * sigG * sigG * C)
            /
            (mu0_vac * sigG * sigG * n_1 * C - (eps2 * n_1 * C_eps + eps1) * ((long double) 1 + n_1 * C));
}



std::ostream & operator << (std::ostream & out, const std::vector<long double> & v) {
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




