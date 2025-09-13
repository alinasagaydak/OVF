#include <iostream>
#include <cmath>
#include <iomanip>

int main() {
    int n_1 = 1;
    while(1.f + powf(2, -n_1) != 1.f) n_1++;
    float eps_1 = powf(2, -(n_1-1));
    int w_1 = 8 * sizeof(float) - (n_1-1) -1;
    int Emin_1 = -pow(2, w_1 - 1) + 2;
    int Emax_1 = pow(2, w_1 - 1) - 1;

    int n_2 = 1;
    while(1. + pow(2, -n_2) != 1.) n_2++;
    double eps_2 = pow(2, -(n_2-1));
    int w_2 = 8 * sizeof(double) - (n_2-1) -1;
    int Emin_2 = -pow(2, w_2 - 1) + 2;
    int Emax_2 = pow(2, w_2 - 1) - 1;
    
    std::cout << "////Single-precision////" << std::endl;
    std::cout << "eps:\t" << eps_1 << std::endl;
    std::cout << "Mantissa bit depth:\t" << n_1-1 << std::endl;
    std::cout << "Minimum binary order:\t" << Emin_1 << std::endl;
    std::cout << "Maximum binary order:\t" << Emax_1 << std::endl;

    std::cout << std::endl;

    std::cout << "////Double-precision////" << std::endl;
    std::cout << "eps:\t" << eps_2 << std::endl;
    std::cout << "Mantissa bit depth:\t" << n_2-1 << std::endl;
    std::cout << "Minimum binary order:\t" << Emin_2 << std::endl;
    std::cout << "Maximum binary order:\t" << Emax_2 << std::endl;

    // Доп вопрос
    std::cout << std::setprecision(64) << 1.f + eps_1 << std::endl;
    std::cout << 1.f + eps_1 + 0.5f*eps_1 << std::endl;
    std::cout << 1.f + 0.5f*eps_1 + eps_1 << std::endl;
    return 0;
}
