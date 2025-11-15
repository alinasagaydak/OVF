#include <iostream>
#include <array>

#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLegend.h"

double func(double x, double u) {
    return -u;
}

// iter functions
double iter_methodEuler(double y_n, double h, double (*func)(double, double), double x_n) {
    double y_n1 = y_n + h * func(x_n, y_n);
    return y_n1;
}

double iter_Runge_Kutta_2(double y_n, double h, double alpha, double (*func)(double, double), double x_n) {
    double y_n1 = y_n + h * ((1-alpha)*func(x_n, y_n) + alpha*func(x_n + h/2./alpha, y_n + h/2./alpha*func(x_n, y_n)));
    return y_n1;
}

double iter_Runge_Kutta_4(double y_n, double h, double (*func)(double, double), double x_n) {
    double k1 = func(x_n, y_n);
    double k2 = func(x_n + h/2., y_n + h/2.*k1);
    double k3 = func(x_n + h/2., y_n + h/2.*k2);
    double k4 = func(x_n + h, y_n + h*k3);
    double y_n1 = y_n + h/6. * (k1 + 2*k2 + 2*k3 + k4);
    return y_n1;
}

// methods
double* methodEuler(const int N, double x_0, double y_0, double h, double (*func)(double, double)) {
    double* arr_y = new double[N + 1];
    for (int i = 0; i < N + 1; i++) {
        arr_y[i] = y_0;
        double y_1 = iter_methodEuler(y_0, h, func, x_0);
        y_0 = y_1;
        x_0 += h;    
    }
    return arr_y;
}

double* Runge_Kutta_2(const int N, double x_0, double y_0, double h, double alpha, double (*func)(double, double)) {
    double* arr_y = new double[N + 1];
    for(int i = 0; i < N + 1; i++) {
        arr_y[i] = y_0;
        double y_1 = iter_Runge_Kutta_2(y_0, h, alpha, func, x_0);
        y_0 = y_1;
        x_0 += h;
    }
    return arr_y;
}

double* Runge_Kutta_4(const int N, double x_0, double y_0, double h, double (*func)(double, double)) {
    double* arr_y = new double[N + 1];
    for(int i = 0; i < N + 1; i++) {
        arr_y[i] = y_0;
        double y_1 = iter_Runge_Kutta_4(y_0, h, func, x_0);
        y_0 = y_1;
        x_0 += h;
    }
    return arr_y;
}


void task6() {
    const int N = 20; // кол-во разбиений
    double x_i = 0.;
    double x_f = 3.;
    double h = (x_f - x_i) / N;
    double alpha = 3./4;

// the Cauchy problem
    double x_0 = 0;
    double u_0 = 1;
   
// exact solution    
    TF1* exactSol = new TF1("exactSol", "TMath::Exp(-x)", 0., 3.);
    exactSol->SetLineColor(kRed);
    exactSol->SetLineWidth(2);

// approximate methods
    double arr_x[N+1];
    for (int i = 0; i < N + 1; i++) arr_x[i] = x_0 + i * h;
    double* arr_y_elr = methodEuler(N, x_0, u_0, h, func);
    double* arr_y_RK2 = Runge_Kutta_2(N, x_0, u_0, h, alpha, func);
    double* arr_y_RK4 = Runge_Kutta_4(N, x_0, u_0, h, func);
    
    TGraph* gr_elr = new TGraph(N+1, arr_x, arr_y_elr);
    TGraph* gr_RK2 = new TGraph(N+1, arr_x, arr_y_RK2);
    TGraph* gr_RK4 = new TGraph(N+1, arr_x, arr_y_RK4);

// The discrepancy
    const int arr_N[3] = {5, 10, 20};
    double discr_elr[3];
    double discr_RK2[3];
    double discr_RK4[3];

    for(int i = 0; i < 3; i++) {
        double tmp_h = (x_f - x_i) / arr_N[i];
        double tmp_x[arr_N[i]+1];
        for (int k = 0; k < arr_N[i] + 1; k++) tmp_x[k] = x_0 + k *tmp_h;
        double* tmp_y_elr = methodEuler(arr_N[i], x_0, u_0, tmp_h, func);
        double* tmp_y_RK2 = Runge_Kutta_2(arr_N[i], x_0, u_0, tmp_h, alpha, func);
        double* tmp_y_RK4 = Runge_Kutta_4(arr_N[i], x_0, u_0, tmp_h, func);

        double max_elr = 0;
        double max_RK2 = 0;
        double max_RK4 = 0;

        for (int j = 0; j < arr_N[i]+1; j++) {
            if (abs(exp(-tmp_x[j]) - tmp_y_elr[j]) > max_elr) max_elr = abs(exp(-tmp_x[j]) - tmp_y_elr[j]);
            if (abs(exp(-tmp_x[j]) - tmp_y_RK2[j]) > max_RK2) max_RK2 = abs(exp(-tmp_x[j]) - tmp_y_RK2[j]);
            if (abs(exp(-tmp_x[j]) - tmp_y_RK4[j]) > max_RK4) max_RK4 = abs(exp(-tmp_x[j]) - tmp_y_RK4[j]);
        }

        discr_elr[i] = max_elr;
        discr_RK2[i] = max_RK2;
        discr_RK4[i] = max_RK4;
    }

    for (auto val : arr_N) std::cout << val << "\t"; std::cout << std::endl;
    std::cout << "Method Euler" << std::endl;
    for (auto val : discr_elr) std::cout << val << "\t"; std::cout << std::endl;
    std::cout << "The Runge-Kutta method of the 2nd order" << std::endl;
    for (auto val : discr_RK2) std::cout << val << "\t"; std::cout << std::endl;
    std::cout << "The Runge-Kutta method of the 4th order" << std::endl;
    for (auto val : discr_RK4) std::cout << val << "\t"; std::cout << std::endl;

// draw results    
    TCanvas* c1 = new TCanvas("c1", "", 1600, 800);
    gStyle->SetLegendFont(4);
    gStyle->SetLegendTextSize(14);
    c1->Divide(3, 1);

// Draw Method Euler    
    c1->cd(1); gPad->SetGrid();
    gr_elr->SetTitle("Method Euler;x;y");
    gr_elr->SetLineWidth(2);
    gr_elr->SetLineColor(kBlue);
    gr_elr->SetMarkerStyle(22);
    gr_elr->Draw("ALP");
    exactSol->Draw("same");
    
    auto lgd1 = new TLegend(0.5, 0.7, 0.8, 0.8);
    lgd1->AddEntry(gr_elr, "Method Euler");
    lgd1->AddEntry(exactSol, "Exact solution");
    lgd1->Draw();

// Draw the Runge-Kutta method of the 2nd order
    c1->cd(2); gPad->SetGrid();
    gr_RK2->SetTitle("The Runge-Kutta method of the 2nd order;x;y");
    gr_RK2->SetLineWidth(2);
    gr_RK2->SetLineColor(kGreen);
    gr_RK2->SetMarkerStyle(22);
    gr_RK2->Draw("ALP");
    exactSol->Draw("same");

    auto lgd2 = new TLegend(0.5, 0.7, 0.82, 0.8);
    lgd2->AddEntry(gr_RK2, "Runge-Kutta 2nd order");
    lgd2->AddEntry(exactSol, "Exact solution");
    lgd2->Draw();

// Draw the Runge-Kutta method of the 4th order
    c1->cd(3); gPad->SetGrid();
    gr_RK4->SetTitle("The Runge-Kutta method of the 4th order;x;y");
    gr_RK4->SetLineWidth(2);
    gr_RK4->SetLineColor(kYellow);
    gr_RK4->SetMarkerStyle(22);
    gr_RK4->Draw("ALP");
    exactSol->Draw("same");
    
    auto lgd3 = new TLegend(0.5, 0.7, 0.81, 0.8);
    lgd3->AddEntry(gr_RK4, "Runge-Kutta 4th order");
    lgd3->AddEntry(exactSol, "Exact solution");
    lgd3->Draw();
    

    delete[] arr_y_elr;
    delete[] arr_y_RK2;
    delete[] arr_y_RK4;
}
