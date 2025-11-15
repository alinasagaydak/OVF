#include <iostream>
#include <math.h>
#include <array>

double func(double x) {
    return 1. / (1. + x*x);
}

double func_erf(double x) {
    return exp(-x*x);
}

double leftRectangle(double a, double b, double N, double (*func)(double)) {
    double h = (b - a) / N;
    double x0 = a;
    double res = 0.;
    for (int i = 0; i < N; i++) {
        res += func(x0) * h;
        x0 += h;
    }
    return res;
}

double rightRectangle(double a, double b, double N, double (*func)(double)) {
    double h = (b - a) / N;
    double x1 = a + h;
    double res = 0;
    for (int i = 0; i < N; i++) {
        res += func(x1) * h;
        x1 += h;
    }
    return res;
}

double trapezoid(double a, double b, double N, double (*func)(double)) {
    double h = (b - a) / N;
    double x0 = a;
    double x1 = a + h;
    double res = 0;
    for (int i = 0; i < N ; i++) {
        res += 0.5 * (func(x0) + func(x1)) * h;
        x0 += h;
        x1 += h;
    }
    return res;
}

double average(double a, double b, double N, double (*func)(double)) {
    double h = (b - a) / N;
    double x0 = a;
    double x1 = a + h;
    double res = 0;
    for (int i = 0; i < N; i++) {
        res += func(0.5*(x0 + x1)) * h;
        x0 += h;
        x1 += h;
    }
    return res;
}

double Simpson(double a, double b, double N, double (*func)(double)) {
    double h = (b - a) / N;
    double x0 = a;
    double x1 = a + h;
    double res = 0;
    for (int i = 0; i < N; i++) {
        res += (h / 6.) * (func(x0) + 4 * func(0.5*(x0 + x1)) + func(x1));
        x0 += h;
        x1 += h;
    }
    return res;
}

double calcErf(double x, double N, double (*Simpson)(double, double, double, double(*func)(double))) {
    double a = 0;
    double b = x;
    double erf = (2. / sqrt(M_PI)) * Simpson(a, b, N, func_erf);
    return erf;
}

double calcErrors(double a, double b, double N, double(*calcInt)(double, double, double, double(*func)(double)), double(*func)(double)) {
    double err = abs(M_PI / 2. - calcInt(a, b, N, func));
    return err;
}

int main() {
    double a = -1.;
    double b = 1.;
    double N = 10.;
    double x = 1.;

    double res_leftRectangle = leftRectangle(a, b, N, func);
    double res_rightRectangle = rightRectangle(a, b, N, func);
    double res_trapezoid = trapezoid(a, b, N, func);
    double res_average = average(a, b, N, func);
    double res_Simpson = Simpson(a, b, N, func);
    double res_erf = calcErf(x, N, Simpson);
    
    std::cout << "left rectangle:\t" << res_leftRectangle << std::endl;
    std::cout << "right rectangle:\t" << res_rightRectangle << std::endl;
    std::cout << "trapezoid:\t" << res_trapezoid << std::endl;
    std::cout << "average:\t" << res_average << std::endl;
    std::cout << "Simpson:\t" << res_Simpson << std::endl;
    std::cout << "erf:\t" << res_erf << std::endl;
    
// Calc Errors
    double Nerr[3] = {10, 20, 40};
    double err_leftRec[3];
    double err_rightRec[3];
    double err_trap[3];
    double err_average[3];
    double err_Simpson[3];
    
    for (int i = 0; i < 3; i++) {
        err_leftRec[i] = calcErrors(a, b, Nerr[i], leftRectangle, func);
        err_rightRec[i] = calcErrors(a, b, Nerr[i], rightRectangle, func);
        err_trap[i] = calcErrors(a, b, Nerr[i], trapezoid, func);
        err_average[i] = calcErrors(a, b, Nerr[i], average, func);
        err_Simpson[i] = calcErrors(a, b, Nerr[i], Simpson, func);
    }    
    
    std::cout << "\n left rectangle:\n" << err_leftRec[0] << "\t" << err_leftRec[1] << "\t" << err_leftRec[2] << std::endl;
    std::cout << "R(10)/R(20): " << err_leftRec[0] / err_leftRec[1] << "\t R(20)/R(40): " << err_leftRec[1] / err_leftRec[2] << std::endl;
    
    std::cout << "\n right rectangle:\n" << err_rightRec[0] << "\t" << err_rightRec[1] << "\t" << err_rightRec[2] << std::endl;
    std::cout << "R(10)/R(20): " << err_rightRec[0] / err_rightRec[1] << "\t R(20)/R(40): " << err_rightRec[1] / err_rightRec[2] << std::endl;
    
    std::cout << "\n trapezoid:\n" << err_trap[0] << "\t" << err_trap[1] << "\t" << err_trap[2] << std::endl;
    std::cout << "R(10)/R(20): " << err_trap[0] / err_trap[1] << "\t R(20)/R(40): " << err_trap[1] / err_trap[2] << std::endl;
    
    std::cout << "\n average:\n" << err_average[0] << "\t" << err_average[1] << "\t" << err_average[2] << std::endl;
    std::cout << "R(10)/R(20): " << err_average[0] / err_average[1] << "\t R(20)/R(40): " << err_average[1] / err_average[2] << std::endl;
    
    std::cout << "\n Simpson:\n" << err_Simpson[0] << "\t" << err_Simpson[1] << "\t" << err_Simpson[2] << std::endl;
    std::cout << "R(10)/R(20): " << err_Simpson[0] / err_Simpson[1] << "\t R(20)/R(40): " << err_Simpson[1] / err_Simpson[2] << std::endl;

    return 0;
}
