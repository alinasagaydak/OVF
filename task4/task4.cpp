#include <iostream>
#include <math.h>
#include <bits/stdc++.h>
#include <iomanip>
#include <functional>

double J0(double x, double t) {
    return 1. / M_PI * cos(x*sin(t));
}

double J1(double x, double t) {
    return 1. / M_PI * cos(t - x*sin(t));
}

double trapezoid(double N, double x, double (*func)(double, double), double a = 0, double b = M_PI, bool deriv_flag = 0) { 
    double h = (b - a) / N; 
    double x0 = a; 
    double x1 = a + h; 
    double res = 0; 
    for (int i = 0; i < N; i++) {
        res += 0.5 * (func(x, x0) + func(x, x1)) * h; 
        x0 += h; 
        x1 += h;
    } 
    return res; 
}

double Simpson(double N, double x, double (*func)(double, double),double a = 0, double b = M_PI) {
    double h = (b - a) / N;
    double x0 = a;
    double x1 = a + h;
    double res = 0;
    for (int i = 0; i < N; i++) {
        res += (h / 6.) * (func(x, x0) + 4 * func(x, 0.5*(x0 + x1)) + func(x, x1));
        x0 += h;
        x1 += h;
    }
    return res;
}


double checkEquality_trap(double N, double x, double (*J0)(double, double), double (*J1)(double, double)) {
    double h = M_PI / N;
    double deriv_J0 =(trapezoid(N, x+h, J0) - trapezoid(N, x-h, J0) ) / 2. / h;
    double int_J1 = trapezoid(N, x, J1);
    std::cout << "deriv_J0 = " << deriv_J0 << std::endl; 
    std::cout << "int_J1 = " << int_J1 << std::endl; 
    return deriv_J0 + int_J1;
}

double checkEquality_simpson(double N, double x, double (*J0)(double, double), double (*J1)(double, double)) {
    double h = M_PI / N;
    double deriv_J0 =(Simpson(N, x+h, J0) - Simpson(N, x-h, J0) ) / 2. / h;
    double int_J1 = Simpson(N, x, J1);
    std::cout << "deriv_J0 = " << deriv_J0 << std::endl; 
    std::cout << "int_J1 = " << int_J1 << std::endl; 
    return deriv_J0 + int_J1;
}

int main() {
    std::srand(time(0));
    double x = (double)rand() / RAND_MAX * 2 * M_PI;
    int precision = 10;
    double N = pow(10, 5);
    std::cout << "x = " << x << std::endl;
    std::cout << "trapezoid: " << checkEquality_trap(N, x, J0, J1) << std::endl;
    std::cout << "Simpson: " << checkEquality_simpson(N, x, J0, J1) << std::endl;

}
