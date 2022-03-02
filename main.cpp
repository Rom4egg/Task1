#include <iostream>
#include "cmath"
#include <vector>
#include "fstream"

constexpr double f = 0;
constexpr double g = 0;

double findZ1(int m, int n, double h, double tao, const std::vector<double> &z1, int pointsH) {
    if (n == 0) return std::sin(m * h) + std::cos(m * h);
    else
        return (f + g - (z1[(n - 1) * pointsH + m] - z1[(n - 1) * pointsH + m - 1]) / h) * tao +
               z1[(n - 1) * pointsH + m];
}

double findZ2(int m, int n, double h, double tao, const std::vector<double> &z2, int pointsH) {
    if (n == 0) return std::sin(m * h) - std::cos(m * h);
    else
        return (f - g + (z2[(n - 1) * pointsH + m + 1] - z2[(n - 1) * pointsH + m]) / h) * tao +
               z2[(n - 1) * pointsH + m];
}

int main() {
    double T = 10;
    int stepsH = 100;
    int pointsH = stepsH + 1;
    int stepsT = 1000;
    int pointsT = stepsT + 1;
    double tao = T / stepsT;
    double h = M_PI / stepsH;

    std::vector<double> z1(pointsH * pointsT, 0);
    std::vector<double> z2(pointsH * pointsT, 0);
    for (int m = 0; m < pointsH; ++m) {
        z1[m] = findZ1(m, 0, h, tao, z1, pointsH);
        z2[m] = findZ2(m, 0, h, tao, z2, pointsH);
    }
    for (int n = 1; n < pointsT; ++n) {
        for (int m = 1; m < pointsH; ++m) {
            z1[n * pointsH + m] = findZ1(m, n, h, tao, z1, pointsH);
        }
        z1[n * pointsH] = std::sin(-n * tao) + std::cos(-n * tao);
    }
    for (int n = 1; n < pointsT; ++n) {
        for (int m = 0; m < pointsH - 1; ++m) {
            z2[n * pointsH + m] = findZ2(m, n, h, tao, z2, pointsH);
        }
        z2[n * pointsH + pointsH - 1] = std::sin(M_PI + n * tao) - cos(M_PI + n * tao);
    }
    std::vector<double> U(pointsH * pointsT);
    std::vector<double> V(pointsH * pointsT);
    for (int i = 0; i < pointsH * pointsT; i++) {
        U[i] = (z1[i] + z2[i]) / 2;
        V[i] = (z1[i] - z2[i]) / 2;
    }
    std::fstream out("gay.csv", std::ios::out);
    out << "t, x, value" << std::endl;
    for (int n = 0; n < pointsT; ++n) {
        for (int m = 0; m < pointsH; ++m) {
            out << n * tao << ", " << m * h << ", " << U[n * pointsH + m] << ", " << std::endl;
        }
    }
    return 0;
}