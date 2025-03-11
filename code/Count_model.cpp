#include "Count_model.h"

Count_Model::Count_Model(long double pos) : z_0(vawe_len* pos) {};

std::complex<long double> U(const Count_Model& model,const long double & z,const long double & z0,const long double & lambda) {
	std::complex<long double> n_1 = sqrt(lambda * lambda - k_1 * k_1);
	std::complex<long double> n_2 = sqrt(lambda * lambda - k_2 * k_2);
	std::complex<long double> n_3 = sqrt(lambda * lambda - k_3 * k_3);

	std::complex<long double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(long double)2 * n_2 * z0);
	std::complex<long double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(long double)2 * n_2 * z0);

	std::complex<long double> Fug = (n_2 * (R - (long double)1) * (r + (long double)1))
		/
		((eps2 * n_1 * (R + (long double)1) - eps1 * n_2 * (R - (long double)1)) * (n_2 * ((long double)1 - r) + n_1 * ((long double)1 + r)) - mu0_vac * sigG * sigG * n_1 * n_2 * (R - (long double)1) * (r + (long double)1))
		*
		(long double)2 * std::complex<long double>(0, 1)
		/
		(omega)
		*
		sigG;

	return Fug * exp(-n_1 * (z + z0)) / n_1;
}

std::complex<long double> U1(const Count_Model& model,const long double & z,const long double & z0,const long double & lambda) {
	std::complex<long double> n_1 = sqrt(lambda * lambda - k_1 * k_1);
	std::complex<long double> n_2 = sqrt(lambda * lambda - k_2 * k_2);
	std::complex<long double> n_3 = sqrt(lambda * lambda - k_3 * k_3);

	std::complex<long double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(long double)2 * n_2 * z0);
	std::complex<long double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(long double)2 * n_2 * z0);

	std::complex<long double> Fug = (n_2 * (R - (long double)1) * (r + (long double)1))
		/
		((eps2 * n_1 * (R + (long double)1) - eps1 * n_2 * (R - (long double)1)) * (n_2 * ((long double)1 - r) + n_1 * ((long double)1 + r)) - mu0_vac * sigG * sigG * n_1 * n_2 * (R - (long double)1) * (r + (long double)1))
		*
		///тут ошибка!
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

	std::complex<long double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(long double)2 * n_2 * z0);
	std::complex<long double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(long double)2 * n_2 * z0);

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

	std::complex<long double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(long double)2 * n_2 * z0);
	std::complex<long double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(long double)2 * n_2 * z0);

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

	std::complex<long double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(long double)2 * n_2 * z0);
	std::complex<long double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(long double)2 * n_2 * z0);

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

	std::complex<long double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(long double)2 * n_2 * z0);
	std::complex<long double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(long double)2 * n_2 * z0);

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

	std::complex<long double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(long double)2 * n_2 * z0);
	std::complex<long double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(long double)2 * n_2 * z0);

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

	std::complex<long double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(long double)2 * n_2 * z0);
	std::complex<long double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(long double)2 * n_2 * z0);

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

	std::complex<long double> R = (n_2 * eps3 - n_3 * eps2) / (n_2 * eps3 + n_3 * eps2) * exp(-(long double)2 * n_2 * z0);
	std::complex<long double> r = (n_2 - n_3) / (n_2 + n_3) * exp(-(long double)2 * n_2 * z0);

	std::complex<long double> Fug = (n_2 * (R - (long double)1) * (r + (long double)1))
		/
		( (eps2 * n_1 * (R + (long double)1) - eps1 * n_2 * (R - (long double)1)) * (n_2 * ((long double)1 - r) + n_1 * ((long double)1 + r)) -
        (mu0_vac * sigG * sigG * n_1 * n_2 * (R - (long double)1) * (r + (long double)1)))
		*
		(long double)2 * std::complex<long double>(0, 1)
		/
		(omega)
		*
		sigG;
    return Fug;
}



std::ostream & operator << (std::ostream & out, const Count_Model& model) {
    out << "k_1: " << k_1 << std::endl;
    out << "k_2: " << k_2 << std::endl;
    out << "k_3: " << k_3 << std::endl;
    return out;
}


