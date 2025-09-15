#include <iostream>
#include <cmath>
#include <vector>

double func(double A, double x) {
    return 1. / tan(sqrt(A*(1.-x))) - sqrt(1./x - 1.);
}

double derivFunc(double A, double x) {
    return 1. / (sin(sqrt(A * (1-x))) * sin(sqrt(A * (1-x)))) * A / 2. / sqrt(A * (1-x)) + 1. / (2 * sqrt(1./x - 1)) * 1. / (x*x);
}

std::vector<double> methodDichotomy(double a, double b, double A, int precision, double (*func)(double, double)) {
    double delta_0 = 0.5 * (a + b);
    double delta = pow(10., -precision);
    double Nnec = abs(log10(delta_0 / delta) / log10(2)); 
    
    double N = 0;
    double f_a = func(A, a);
    double a_in = a;
    double b_in = b;
    while(fabs(a - b) > pow(10, -precision)) {
        double f_mid = func(A, 0.5 * (a + b));
        if (f_a * f_mid < 0) b = 0.5 * (a + b);
        else a = 0.5 * (a + b);
        
        if (N > 10 * Nnec) {break; return std::vector<double> {0,0};}
        N++;
    }
    std::vector<double> res1 = {0.5 * (a + b), N};
    
    if (fabs(a_in - res1[0]) > fabs(b_in - res1[0])) {
        N = 0;
        //b = res1[0]-0.01;
        b = 0.5 * (a_in + b_in);
        f_a = func(A, a_in);
        while(fabs(a - b) > pow(10, -precision)) {
            double f_mid = func(A, 0.5 * (a + b));
            if (f_a * f_mid < 0) b = 0.5 * (a + b);
            else a = 0.5 * (a + b);
        
            if (N > 10 * Nnec) {break; return std::vector<double> {0,0};}
            N++;
        }
        std::vector<double> res2 = {0.5 * (a + b), N};
        return res1[0] > res2[0] ? res2 : res1;
    }
    return res1;
}

std::vector<double> methodSimpleIter(double A, double a, double b, double x_0, int precision, double (*func)(double, double), double (*derivFunc)(double, double)) {
    double N = 0;
    double x_1 = a;
    double a_in = a;
    double b_in = b;
    double lambda = 1. / derivFunc(A, x_0);
    double q = 1 - lambda * derivFunc(A, x_0);
    while(fabs((x_1 - x_0) / (1.q - 1)) > pow(10, -precision)) {
        double tmp = x_1;
        double sgn = derivFunc(A, x_0) > 0 ? 1. : -1.;
        x_1 = x_0 - sgn * lambda * func(A, x_0);
        x_0 = tmp;
        N++;
        if (N > pow(10,4)) {break; return std::vector<double> {0,0};}
        
        sgn = derivFunc(A, x_0) > 0 ? 1. : -1.;
        q = 1 - sgn * lambda * derivFunc(A, x_0);
    }
    std::vector<double> res1 = {0.5*(x_0 + x_1), N};

    if (fabs(a_in - res1[0]) > fabs(b_in - res1[0])) {
        N = 0;
        x_1 = a_in;
        //x_0 = res1[0] - 0.1;
        x_0 = 0.5 * (a_in + b_in);
        q = 1 - lambda * derivFunc(A, x_0);
        while(fabs((x_1 - x_0) / (1./q - 1)) > pow(10, -precision)) {
            double tmp = x_1;
            double sgn = derivFunc(A, x_0) > 0 ? 1. : -1.;
            x_1 = x_0 - sgn * lambda * func(A, x_0);
            x_0 = tmp;
            N++;
            if (N > 10000) {break; return std::vector<double> {0,0};}

            sgn = derivFunc(A, x_0) > 0 ? 1. : -1.;
            q = 1 - sgn * lambda * derivFunc(A, x_0);
        }
        std::vector<double> res2 = {0.5*(x_0 + x_1), N};
        return res1[0] > res2[0] ? res2 : res1;
    }
    return res1;
}

std::vector<double> methodNewton(double A, double a, double b, double x_0, int precision, double (*func)(double, double), double (*derivFunc)(double, double)) {
    double N = 0;
    double x_1 = a;
    double a_in = a;
    double b_in = b;
    while(fabs(x_1 - x_0) > pow(10, -precision)) {
        double tmp = x_1;
        x_1 = x_0 - func(A, x_0) / derivFunc(A, x_0);
        x_0 = tmp;
        N++;
        if (N > 10000) {break; return std::vector<double> {0,0};}
    }
    std::vector<double> res1 = {0.5 * (x_0 + x_1), N};
    
    if (fabs(a_in - res1[0]) > fabs(b_in - res1[0])) {
        N = 0;
        x_1 = a_in;
        //x_0 = res1[0] - 0.1;
        x_0 = 0.5 * (a_in + b_in);
        while(fabs(x_1 - x_0) > pow(10, -precision)) {
            double tmp = x_1;
            x_1 = x_0 - func(A, x_0) / derivFunc(A, x_0);
            x_0 = tmp;
            N++;
            if (N > 10000) {break; return std::vector<double> {0,0};}
        }
        std::vector<double> res2 = {0.5 * (x_0 + x_1), N};
        return res1[0] > res2[0] ? res2 : res1;
    }
    
    return res1;
}

int main() {
    double a = 0.01;
    double b = 0.99;
    double A = 25.;
    int precision = 5;

    std::vector<double> xD = methodDichotomy(a, b, A, precision, func);
    std::vector<double> xI = methodSimpleIter(A, a, b, b, precision, func, derivFunc);
    std::vector<double> xN = methodNewton(A, a, b, b, precision, func, derivFunc);

    std::cout << "The dichotomy method: " << xD[0] << "\tN: " << xD[1] << std::endl;
    std::cout << "The method of simple iterations: " << xI[0] << "\tN: " << xI[1] << std::endl;
    std::cout << "Newton's method: " << xN[0] << "\tN: " << xN[1] << std::endl;
    return 0;
}
