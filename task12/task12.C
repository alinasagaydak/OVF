#include <TMath.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <iostream>
#include <string>
#include <set>
#include <iomanip>
#include <cmath>

#include <stdio.h>
#include <math.h>

void direct_FFT(int N, double* signal, double* transformed_real, double* transformed_imag) {
    for (int i = 0; i < N; i++) {
        double real_tmp = 0.0;
        double imag_tmp = 0.0;
        for (int k = 0; k < N; k++) {
            double angle = 2*M_PI*i*k/N;
            real_tmp += signal[k] * cos(angle);
            imag_tmp += signal[k] * sin(angle);
        }
        transformed_real[i] = real_tmp / N;
        transformed_imag[i] = imag_tmp / N;
    }
}

void reverse_FFT(int N, double* transformed_real, double* transformed_imag, double* restored_real) {
    for (int i = 0; i < N; i++) {
        double real_tmp = 0.0;
        for (int k = 0; k < N; k++) {
            double angle = -2*M_PI*i*k/N;
            real_tmp += transformed_real[k] * cos(angle) - transformed_imag[k] * sin(angle);
        }
        restored_real[i] = real_tmp;
    }

}

void task12() {
    const int N = 1000;
    double a0 = 1.0;
    double a1 = 0.002;
    double w0 = 5.1;
    double w1 = 25.5;
    double t0 = 0.0;
    double t1 = 2*M_PI;
    
    double t[N];
    double signal[N];
    
    double transformed_real[N];  
    double transformed_imag[N];  
    double restored_real[N];    

    double transformed_hann_real[N];
    double transformed_hann_imag[N];
    double restored_hann_real[N];
    
    double hann_window[N];
    
    for (int i = 0; i < N; i++) {
        t[i] = t0 + i * (t1 - t0) / (N - 1);
        signal[i] = a0 * sin(w0 * t[i]) + a1 * sin(w1 * t[i]);
    }

    direct_FFT(N, signal, transformed_real, transformed_imag);
    reverse_FFT(N, transformed_real, transformed_imag, restored_real);
    
    // Hann
    for (int k = 0; k < N; k++) {
        hann_window[k] = 0.5 * (1 - cos(2 * M_PI * k / N));
    }

    double hann_signal[N];
    for (int i = 0; i < N; i++) hann_signal[i] = signal[i] * hann_window[i];
   
    direct_FFT(N, hann_signal, transformed_hann_real, transformed_hann_imag);
    reverse_FFT(N, transformed_hann_real, transformed_hann_imag, restored_hann_real);
    
   
    // Power
    double power_no_window[N];
    double power_hann[N];
    double frequency[N];
    for (int i = 0; i < N; i++) {
        power_no_window[i] = transformed_real[i] * transformed_real[i] + 
                                transformed_imag[i] * transformed_imag[i];
        power_hann[i] = transformed_hann_real[i] * transformed_hann_real[i] + 
                           transformed_hann_imag[i] * transformed_hann_imag[i];
        frequency[i]= i * (1./(t[1]-t[0]))/N *2*M_PI;
    }

    double log_power_no_window[N];
    double log_power_hann[N];
    for (int i = 0; i < N; i++) {
        log_power_no_window[i] = log(power_no_window[i]);
        log_power_hann[i] = log(power_hann[i]);
    }
    

    TCanvas *canvas = new TCanvas();
    canvas->Divide(3,1);

    canvas->cd(1);
    TGraph* real_signal = new TGraph(N, t, signal);
    real_signal->SetLineWidth(2);
    real_signal->SetTitle("real signal;t;f(t)");
    real_signal->Draw();
    canvas->cd(2);
    TGraph* restored_signal = new TGraph(N, t, restored_real);
    restored_signal->SetLineWidth(2);
    restored_signal->SetTitle("restored signal;t;f(t)");
    restored_signal->Draw();
    canvas->cd(3);
    TGraph* restored_hann = new TGraph(N, t, restored_hann_real);
    restored_hann->SetLineWidth(2);
    restored_hann->SetTitle("restored hann signal;t;f(t)");
    restored_hann->Draw();


    TCanvas *canvas1 = new TCanvas();
    canvas1->Divide(2,1);

    canvas1->cd(1);
    TGraph* power_1 = new TGraph(N/2, frequency, log_power_no_window);
    power_1->SetLineWidth(2);
    power_1->SetTitle("no window;frequency;power");
    power_1->Draw();
    canvas1->cd(2);
    TGraph* power_2 = new TGraph(N/2, frequency, log_power_hann);
    power_2->SetLineWidth(2);
    power_2->SetTitle("hann window;frequency;power");
    power_2->Draw();
}
