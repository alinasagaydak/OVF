#include <iostream>
#include <array>
#include <vector>

double func1(double* coeffs, double t, double x, double y) {
    return coeffs[0] * x - coeffs[1] * x * y;
}

double func2(double* coeffs, double t, double x, double y) {
    return coeffs[2] * x * y - coeffs[3] * y;
}

std::vector<double> iter_Runge_Kutta_2
(double* coeffs, double t_n, double x_n, double y_n, double h, double alpha, double (*func1)(double*, double, double, double), double (*func2)(double*, double, double, double)) {
    double delta_t = h/2./alpha;
    double delta_x = h/2./alpha * func1(coeffs, t_n, x_n, y_n);
    double delta_y = h/2./alpha * func2(coeffs, t_n, x_n, y_n);
    double x_n1 = x_n + h * ((1-alpha) * func1(coeffs, t_n, x_n, y_n) + alpha * func1(coeffs, t_n + delta_t, x_n + delta_x, y_n + delta_y));
    double y_n1 = y_n + h * ((1-alpha) * func2(coeffs, t_n, x_n, y_n) + alpha * func2(coeffs, t_n + delta_t, x_n + delta_x, y_n + delta_y));
    std::vector<double> res = {x_n1, y_n1};
    return res;
}

std::vector<std::vector<double>> Runge_Kutta_2
(const int N, double* coeffs, double t_0, double x_0, double y_0, double h, double alpha, double (*func1)(double*, double, double, double), double (*func2)(double*, double, double, double)) {
    std::vector<std::vector<double>> vec_res;
    for(int i = 0; i < N+1; i++) {
        vec_res.push_back({x_0, y_0});
        double x_1 = iter_Runge_Kutta_2(coeffs, t_0, x_0, y_0, h, alpha, func1, func2)[0];
        double y_1 = iter_Runge_Kutta_2(coeffs, t_0, x_0, y_0, h, alpha, func1, func2)[1];
        x_0 = x_1;
        y_0 = y_1;
        t_0 += h;
    }
    return vec_res;
}

void task7() {
    double coeffs[4] = {10, 2, 2, 10}; // a b c d
    const int N = 150;
    double t_i = 0.;
    double t_f = 10.;
    double h = (t_f - t_i) / N;
    double alpha = 3./4;

    double t_0 = 0.;
    double x_0 = 2; // жертвы
    double y_0 = 2; // хищники

    double arr_t[N+1];
    for (int i = 0; i < N + 1; i++) arr_t[i] = t_0 + h * i;

    std::vector<std::vector<double>> vec_xy = Runge_Kutta_2(N, coeffs, t_0, x_0, y_0, h, alpha, func1, func2);
    double arr_x[N+1];
    double arr_y[N+1];
    for (int i = 0; i < vec_xy.size(); i++) {
        arr_x[i] = vec_xy[i][0];
        arr_y[i] = vec_xy[i][1];
    }

    TCanvas* c1 = new TCanvas("c1", "", 1600, 800);
    c1->Divide(2, 1);

    c1->cd(1); gPad->SetGrid();
    TGraph* phaseTr = new TGraph(N+1, arr_x, arr_y);
    phaseTr->SetTitle("phase trajectory; x(t);y(t)");
    phaseTr->SetLineColor(kBlue);
    phaseTr->Draw("AC");

    c1->cd(2); gPad->SetGrid();
    TGraph* x_t = new TGraph(N+1, arr_t, arr_x);
    TGraph* y_t = new TGraph(N+1, arr_t, arr_y);
    x_t->SetLineColor(kRed); x_t->SetLineWidth(2);
    y_t->SetLineColor(kBlack); y_t->SetLineWidth(2);
    TMultiGraph* mgr = new TMultiGraph();
    mgr->Add(x_t);
    mgr->Add(y_t);
    mgr->SetTitle("x(t) and y(t);t;x(t) / y(t)");
    mgr->Draw("AC");

    auto lgd = new TLegend(0.75, 0.75, 0.85, 0.85);
    lgd->AddEntry(x_t, "x(t)");
    lgd->AddEntry(y_t, "y(t)");
    lgd->Draw();
}
