#include <iostream>
#include <math.h>
#include <iomanip>

double funcSum_StL(int n) {
    double sum = 0;
    for (int i = 1; i < n+1; i++) {
        sum += pow((-1.),i) / i;
    }
    return sum;
}

double funcSum_LtS(int n) {
    double sum = 0;
    for (int i = n; i > 0; i--) {
        sum += pow((-1.),i) / i;
    }
    return sum;
}

double funcSum_negative(int n, bool flag) {
    double sum = 0;
    if (flag) { // from small to large
        for (int i = 1; i < n+1; i++) {
            sum += pow((-1.), i) / i;
            i++;    
        }
        return sum;
    }
    if (!flag) { // from large to small
        for (int i = n-1; i > 0; i--) {
            sum += pow((-1.),i) / i;
            i--;
        }
        return sum;
    }
    return -1;
}

double funcSum_positive(int n, bool flag) {
    double sum = 0;
    if (flag) { // from small to large
        for (int i = 2; i < n+1; i++) {
            sum += pow((-1.),i) / i;
            i++;
        }
        return sum;
    }
    if (!flag) { // from large to small
        for (int i = n; i > 0; i--) {
            sum += pow((-1.),i) / i;
            i--;
        }
        return sum;
    }
    return -1;
}

int main() {
    int n = pow(10, 5);
    double sum_small_to_large = funcSum_StL(n);
    double sum_large_to_small = funcSum_LtS(n);
    double sum_posNeg_small_to_large = funcSum_negative(n, 1) + funcSum_positive(n, 1);
    double sum_posNeg_large_to_small = funcSum_negative(n, 0) + funcSum_positive(n, 0);

    std::cout << "Sum from small to large:\t" << std::setprecision(64) << sum_small_to_large << std::endl;
    std::cout << "Sum from large to small:\t" << sum_large_to_small << std::endl;
    std::cout << "Sum posNeg from small to large:\t" << sum_posNeg_small_to_large << std::endl;
    std::cout << "Sum posNeg from large to small:\t" << sum_posNeg_large_to_small << std::endl;

    return 0;
}
