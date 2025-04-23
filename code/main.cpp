#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include "Count_model.h"
#include <fstream>
#include <iostream>
#include <fstream>
#include <thread>
#include <vector>
#include <mutex>
#include <iomanip>
#include <cstdio>
#include <string>



std::complex<double> W(const Count_Model& model, const double & z, const double & z0, const double & lambda);
std::complex<double> U(const Count_Model& model, const double & z, const double & z0, const double & lambda);
std::complex<double> E_rho(const Count_Model&, const double &, const double &, const double &);
std::complex<double> E_phi(const Count_Model&, const double &, const double &, const double &);
std::complex<double> E_z(const Count_Model&, const double &, const double &, const double &);

std::complex<double> H_rho(const Count_Model&, const double &, const double &, const double &);
std::complex<double> H_phi(const Count_Model&, const double &, const double &, const double &);
std::complex<double> H_z(const Count_Model&, const double &, const double &, const double &);
std::complex<double> F_wg(const Count_Model&, const double &, const double &);
std::complex<double> F_ug(const Count_Model&, const double &, const double &);
std::complex<double> F_ug_n_2_is_zero(const Count_Model&, const double &, const double &);
std::complex<double> F_ug_n_1_is_zero(const Count_Model&, const double &, const double &);
std::complex<double> F_wg_n_2_is_zero(const Count_Model&, const double &, const double &);
std::complex<double> F_wg_n_1_is_zero(const Count_Model&, const double &, const double &);
int n = 0;

// Мьютекс для синхронизации доступа к файлу
std::mutex file_mutex;


std::ostream & operator << (std::ostream & out, const std::complex<double> & z) {
    out << "{" << z.real() << ", " << z.imag() << "}";
    return out;
}



void run_and_write_to_file(const std::string& filename,
                            std::complex<double>(*func)(const Count_Model&, const double &, const double &, const double &),
                            double z, double rho,const Count_Model & model, std::string name) {
    std::complex<double> res = func(model, z, model.z_0, rho);
    std::lock_guard<std::mutex> lock(file_mutex); // Защита мьютексом
    n++;
    std::ofstream file(filename, std::ios::app);  // Открытие файла в режиме добавления
    if (file.is_open()) {
        file << std::fixed << std::setprecision(40);
        file << "z: " << z << std::endl;
        file << "R: " << rho << std::endl;
        file << name << std::endl << res << std::endl;
        std::cout << "I have counted " << n << " of 480 components!" << std::endl;
    } else {
        std::cerr << "Не удалось открыть файл: " << filename << std::endl;
    }
}



void count_spot(const Count_Model & model, double z, double rho) {
    std::vector<std::thread> threads;
    threads.emplace_back(run_and_write_to_file, "_stdout.txt", E_rho, z, rho, model, "E_rho");
    threads.emplace_back(run_and_write_to_file, "_stdout.txt", H_rho, z, rho, model, "H_rho");
    threads.emplace_back(run_and_write_to_file, "_stdout.txt", E_phi, z, rho, model, "E_phi");
    threads.emplace_back(run_and_write_to_file, "_stdout.txt", H_phi, z, rho, model, "H_phi");
    threads.emplace_back(run_and_write_to_file, "_stdout.txt", E_z, z, rho, model, "E_z");
    threads.emplace_back(run_and_write_to_file, "_stdout.txt", H_z, z, rho, model, "H_z");


    for (auto& t : threads) {
        if (t.joinable()) {
            t.join();
        }
    }
}


void count_one_comp(const Count_Model & model, std::vector<double> & spots_z, double rho, std::string comp_name,
                    std::complex<double>(*func)(const Count_Model&, const double &, const double &, const double &)) {
    const std::string f_name= comp_name + ".out";
    std::remove(f_name.c_str());
    std::vector<std::thread> threads;
    for (auto z: spots_z) {
        threads.emplace_back(run_and_write_to_file, comp_name + ".out", func, z, rho, model, comp_name);
    }

    for (auto& t : threads) {
        if (t.joinable()) {
            t.join();
        }
    }
}



void count_one_comp_non_paralel(const Count_Model & model, std::vector<double> & spots_z, double rho, std::string comp_name,
                    std::complex<double>(*func)(const Count_Model&, const double &, const double &, const double &)) {
    const std::string f_name= comp_name + ".out";
    std::remove(f_name.c_str());

    for (auto z: spots_z) {
        run_and_write_to_file(comp_name + ".out", func, z, rho, model, comp_name);
    }


}

void count_f_near_spot(const Count_Model & model, std::complex<double>(*func) (const Count_Model&, const double &, const double &),
                       double spot, int n_spots = 100,
                       double step = 0.000001, std::string out_file_name = "f_show_file") {

    std::vector<double> xs(n_spots * 2);
    std::vector<std::complex<double>> f_values(n_spots * 2);

    for (int i = -n_spots; i < n_spots; ++i) {
        /*if (i == 0) {
            xs[i + n_spots] = -1;
            f_values[i + n_spots] = {-1, -1};
            continue;
        }*/
        xs[n_spots + i] = spot + (i * step);
        f_values[n_spots + i] = func(model, xs[n_spots + i], model.z_0);
    }

    std::ofstream out(out_file_name + ".out");
    out << "spots" << std::endl;
    for (auto el : xs) {
        ///out << std::fixed << std::setprecision(40);
        out << el << ' ';
    }
    out << std::endl << "f values real" << std::endl;
    for (auto el : f_values) {
        ///out << std::fixed << std::setprecision(40);
        out << el.real() << ' ';
    }

    out << std::endl << "f values im" << std::endl;
    for (auto el : f_values) {
        ///out << std::fixed << std::setprecision(40);
        out << el.imag() << ' ';
    }

    out << std::endl <<  n_spots << std::endl;

    out.close();

}


void count_f_near_spot_4_args(const Count_Model & model,
                       std::complex<double>(*func) (const Count_Model&, const double &, const double &, const double &),
                       double spot, int n_spots = 100,
                       double step = 0.000001, std::string out_file_name = "f_show_file") {
    std::vector<double> xs(n_spots * 2);
    std::vector<std::complex<double>> f_values(n_spots * 2);
    for (int i = -n_spots; i < n_spots; ++i) {
        if (i == 0) {
            xs[i + n_spots] = -1;
            f_values[i + n_spots] = {-1, -1};
            continue;
        }
        xs[n_spots + i] = spot + (i * step);
        f_values[n_spots + i] = func(model, 0.000001, model.z_0, xs[n_spots + i]);
    }

    std::ofstream out(out_file_name + ".out");

    out << "spots" << std::endl;
    for (auto el : xs) {
        out << std::fixed << std::setprecision(40) << el << ' ';
    }
    out << std::endl << "f values real" << std::endl;
    for (auto el : f_values) {
        out << std::fixed << std::setprecision(40) << el.real() << ' ';
    }

    out << std::endl << "f values im" << std::endl;
    for (auto el : f_values) {
        out << std::fixed << std::setprecision(40) << el.imag() << ' ';
    }

    out << std::endl <<  n_spots << std::endl;

    out.close();

}


std::vector<double> prepare_z(const Count_Model & model, int n = 4) {
    std::vector<double> res;
    /*for (int i = 0; i < n / 4; ++i) {
        res.push_back(0.0 + model.z_0 / ((double) (n / 4)) * i);
    }
    for(int i = 0; i < n / 4; ++i) {
        res.push_back(model.z_0 + model.z_0 / ((double) (n / 4)) * i);
    }
    for(int i = 0; i < n / 4; ++i) {
        res.push_back(model.z_0 + model.z_0 * i);
    }*/

    double k = 10.;
    for(int i = 0; i < n; ++i) {
        for (int j =  1; j < 100; ++j)
            res.push_back(model.z_0 + model.z_0 * k * j / 10.);
        k *= (double) 10.;
    }

    

    return res;
}







int main()
{
    void* staus = freopen("stderr.txt", "w", stderr);
    ///freopen("stdout.txt", "w", stdout);
    Count_Model model(0.2, 0.1, 10000, 10000);

    ///std::cout << model.z_0 << std::endl << prepare_z(model) << std::endl;
    std::vector<double> vec_of_spots = prepare_z(model, 9);
    std::cout << model << '\n';
    count_one_comp_non_paralel(model, vec_of_spots, 0., "E_phi", E_phi);
    /*std::cout << vec_of_spots << std::endl;

    for (auto el : vec_of_spots) {
        try {
            count_spot(model, el, 0.);
        }
        catch (...) {
            std::cerr << "Ошибка в точке z = " << el << " r = 0" << std::endl;
        }
    }*/


    /*for (auto el : vec_of_spots) {
        try {
            count_spot(model, model.z_0, el);
        }
        catch (...) {
            std::cerr << "Ошибка в точке r = " << el << " r = " << model.z_0 << std::endl;
        }
    }*/
    std::cout << model << std::endl;
    /*count_one_comp(model, vec_of_spots, 0., "G_H_z", H_z);
    count_one_comp(model, vec_of_spots, 0., "G_E_z", E_z);
    count_one_comp(model, vec_of_spots, 0., "G_H_rho", H_rho);
    count_one_comp(model, vec_of_spots, 0., "G_H_phi", H_phi);
    count_one_comp(model, vec_of_spots, 0., "G_E_phi", E_phi);
    count_one_comp(model, vec_of_spots, 0., "G_E_rho", E_rho);

    std::cout << "Все функции завершили работу, данные записаны в output.txt.\n";
    */

    //count_f_near_spot(model, F_ug, k_2.real(), 100, 0.00000000000001, "F_ug_n_2_file");
    ///count_f_near_spot(model, F_ug_n_2_is_zero, k_2.real(), 100, 0.0001, "F_ug_n_2_0_file");

    ///std::cout << F_wg(model, k_1.real(), model.z_0) << std::endl;
    ///std::cout << F_wg_n_1_is_zero(model, k_1.real(), model.z_0) << std::endl;

    return 0;

}

