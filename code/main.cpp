// Graphene model 2.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include "Count_model.h"
#include <fstream>
#include <iostream>
#include <fstream>
#include <thread>
#include <vector>
#include <mutex>
#include <iomanip>
std::complex<long double> W(Count_Model& model, long double z, long double z0, long double lambda);
std::complex<long double> U(Count_Model& model, long double z, long double z0, long double lambda);
std::complex<long double> E_rho(const Count_Model&, const long double &, const long double &, const long double &);
std::complex<long double> E_phi(const Count_Model&, const long double &, const long double &, const long double &);
std::complex<long double> E_z(const Count_Model&, const long double &, const long double &, const long double &);

std::complex<long double> H_rho(const Count_Model&, const long double &, const long double &, const long double &);
std::complex<long double> H_phi(const Count_Model&, const long double &, const long double &, const long double &);
std::complex<long double> H_z(const Count_Model&, const long double &, const long double &, const long double &);
std::complex<long double> F_wg(const Count_Model&, const long double &, const long double &);
std::complex<long double> F_ug(const Count_Model&, const long double &, const long double &);
int n = 0;

// Мьютекс для синхронизации доступа к файлу
std::mutex file_mutex;


void run_and_write_to_file(const std::string& filename,
                            std::complex<long double>(*func)(const Count_Model&, const long double &, const long double &, const long double &),
                            long double z, long double rho,const Count_Model & model, std::string name) {
    std::complex<long double> res = func(model, rho, z, model.z_0);
    std::lock_guard<std::mutex> lock(file_mutex); // Защита мьютексом
    n++;
    std::ofstream file(filename, std::ios::app);  // Открытие файла в режиме добавления
    if (file.is_open()) {
        file << "z: " << z << std::endl;
        file << "R: " << rho << std::endl;
        file << name << std::endl << res << std::endl;
        std::cout << "I have counted " << n << " of 480 components!" << std::endl;
    } else {
        std::cerr << "Не удалось открыть файл: " << filename << std::endl;
    }
}



void count_spot(const Count_Model & model, long double z, long double rho) {
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


void count_one_comp(const Count_Model & model, std::vector<long double> & spots_z, long double rho, std::string comp_name,
                    std::complex<long double>(*func)(const Count_Model&, const long double &, const long double &, const long double &)) {
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


void count_f_near_spot(const Count_Model & model, std::complex<long double>(*func) (const Count_Model&, const long double &, const long double &),
                       long double spot, int n_spots = 100,
                       long double step = 0.000001, std::string out_file_name = "f_show_file") {
    std::vector<long double> xs(n_spots * 2);
    std::vector<std::complex<long double>> f_values(n_spots * 2);
    for (int i = -n_spots; i < n_spots; ++i) {
        if (i == 0) continue;
        xs[n_spots + i] = spot + (i * step);
        f_values[n_spots + i] = func(model, xs[n_spots + i], model.z_0);
    }

    std::ofstream out(out_file_name + ".out");

    out << "spots" << std::endl;
    for (auto el : xs) {
        out << std::fixed << std::setprecision(20) << el << ' ';
    }
    out << "f values real" << std::endl;
    for (auto el : f_values) {
        out << std::fixed << std::setprecision(20) << el.real() << ' ';
    }

    out << "f values im" << std::endl;
    for (auto el : f_values) {
        out << std::fixed << std::setprecision(20) << el.imag() << ' ';
    }

    out.close();

}


std::vector<long double> prepare_z(const Count_Model & model, int n = 4) {
    std::vector<long double> res;
    for(int i = 0; i < n / 4; ++i) {
        res.push_back(0.0 + model.z_0 / ((long double) (n / 4)) * i);
    }
    for(int i = 0; i < n / 4; ++i) {
        res.push_back(model.z_0 + model.z_0 / ((long double) (n / 4)) * i);
    }
    for(int i = 0; i < n / 4; ++i) {
        res.push_back(model.z_0 + model.z_0 * i);
    }

    long double k = 10.;
    for(int i = 0; i < n / 4; ++i) {
        res.push_back(model.z_0 + model.z_0 * k);
        k *= (long double) 10.;
    }

    return res;
}



std::ostream & operator << (std::ostream & out, const std::vector<long double> & v) {
    for (auto el : v) {
        out << el << ' ';
    }

    return out;
}



int main()
{
    freopen("stderr.txt", "w", stderr);
    ///freopen("stdout.txt", "w", stdout);
    Count_Model model(0.2);

    ///std::cout << model.z_0 << std::endl << prepare_z(model) << std::endl;
    //std::vector<long double> vec_of_spots = prepare_z(model, 40);
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

    count_f_near_spot(model, F_ug, k_1.real(), 100, 0.00000001, "F_ug_file");

    return 0;

}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
