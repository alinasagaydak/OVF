#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include <cmath>

const int N = 1000;
const double x_min = -10.0;
const double x_max = 10.0;
const int max_iter = 1000;
const double accur = 1e-12;

double potential(double x) {
    return 0.5 * x * x;  
}

double analytic_wave_func(double x) {
    return pow(M_PI, -0.25) * exp(-x * x / 2.0);
}

void matrix_H(double *diag, double *sub, double *sup, double h) {
    for (int i = 0; i < N; i++) {
        double x = x_min + i * h;
        diag[i] = 1./(h*h) + potential(x);
    }
    for (int i = 0; i < N-1; i++) {
        sub[i] = -0.5/(h*h); // lower
        sup[i] = -0.5/(h*h); // upper
    }
}

void solve(const double *diag, const double *sub, const double *sup, const double *b, double *x) {
    double c_prime[N-1]; // x[i] = d_prime[i] - c_prime[i] * x[i+1]
    double d_prime[N];
    c_prime[0] = sup[0] / diag[0];
    d_prime[0] = b[0] / diag[0];
    for (int i = 1; i < N-1; i++) {
        double denom = diag[i] - sub[i-1] * c_prime[i-1];
        c_prime[i] = sup[i] / denom;
        d_prime[i] = (b[i] - sub[i-1] * d_prime[i-1]) / denom;
    }
    d_prime[N-1] = (b[N-1] - sub[N-2] * d_prime[N-2]) / (diag[N-1] - sub[N-2] * c_prime[N-2]);
    x[N-1] = d_prime[N-1];
    for (int i = N-2; i >= 0; i--) {
        x[i] = d_prime[i] - c_prime[i] * x[i+1];
    }
}

double dot_product(const double *a, const double *b, double h) {
    double result = 0.0;
    for (int i = 0; i < N; i++) {
        result += a[i] * b[i];
    }
    return result * h;
}

void normalize(double *vec, double h) {
    double norm_sq = dot_product( vec, vec, h);
    double norm = sqrt(norm_sq);
    for (int i = 0; i < N; i++) {
        vec[i] /= norm;
    }
}

void matrix_vector_mult(const double *diag, const double *sub, const double *sup, const double *vec, double *result) {
    result[0] = diag[0] * vec[0] + sup[0] * vec[1];
    for (int i = 1; i < N-1; i++) {
        result[i] = sub[i-1] * vec[i-1] + diag[i] * vec[i] + sup[i] * vec[i+1];
    }
    result[N-1] = sub[N-2] * vec[N-2] + diag[N-1] * vec[N-1];
}


double inverse_iteration(const double *diag, const double *sub, const double *sup, double *eigenvector, double shift, double h) {
    double diag_shift[N];//сдвиг
    for (int i = 0; i < N; i++) {
        diag_shift[i] = diag[i] - shift;
    }
    
    for (int i = 0; i < N; i++) { //приближение
        eigenvector[i] = 1.;
    }
    normalize(eigenvector, h);
    
    double eigenvalue_old = 0.0;
    double eigenvalue_new;
    
    for (int iter = 0; iter < max_iter; iter++) {
        double psi_new[N];
        solve(diag_shift, sub, sup, eigenvector, psi_new);
        
        memcpy(eigenvector, psi_new, N * sizeof(double));
        normalize(eigenvector, h);
        
        double H_psi[N];
        matrix_vector_mult(diag, sub, sup, eigenvector, H_psi);
        eigenvalue_new = dot_product(eigenvector, H_psi, h);
        
        if (fabs(eigenvalue_new - eigenvalue_old) < accur) break;
        
        eigenvalue_old = eigenvalue_new;
    }
    
    return eigenvalue_new;
}

void task11() {
    double h = (x_max - x_min) / (N - 1);  
    
    double diag[N];     
    double sub[N-1];    
    double sup[N-1];   
    matrix_H(diag, sub, sup, h);
    
    double psi[N];
    
    double shift = 0.4;
    double energy = inverse_iteration(diag, sub, sup, psi, shift, h);

    double psi_analytic[N];
    double x_arr[N];
    for (int i = 0; i < N; i++) {
        double x = x_min + i * h;
        psi_analytic[i] = analytic_wave_func(x);
        x_arr[i] = x;
    }
    
    printf("Численно: %.10f\n", energy);
    printf("Аналитически: 0.5\n");
    printf("Ошибка: %.2e\n", fabs(energy - 0.5));
    
    TCanvas *canvas = new TCanvas();
    TGraph* g = new TGraph(N, x_arr, psi);
    g->SetLineWidth(6);
    g->SetLineColor(kBlue);
    g->Draw();
    TGraph* g1 = new TGraph(N, x_arr, psi_analytic);
    g1->SetLineWidth(2);
    g1->SetLineColor(kRed);
    g1->Draw("same");

    auto lgd = new TLegend(0.6, 0.6, 0.8, 0.7);
    lgd->AddEntry(g, "Numerical solution");
    lgd->AddEntry(g1, "Analytical solution");
    lgd->Draw();
}
