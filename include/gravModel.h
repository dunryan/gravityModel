#include <vector>
#include <cmath>
#include "./constants.h"
#include "../../eigen/Eigen/Dense"

void provideHarmonicCoefficientsJGM3(std::vector<std::vector<double>>& C, std::vector<std::vector<double>>& S);
void provideHarmonicCoefficientsEGM2008(std::vector<std::vector<double>>& C, std::vector<std::vector<double>>& S);
double legendreP(int l, int m, double x);
double dLegendreP(int l, int m, double x);
std::vector<double> computeGravitationalAcceleration(double r, double theta, double lambda, const std::vector<std::vector<double>>& C, const std::vector<std::vector<double>>& S);

std::vector<std::vector<double>> C(1, std::vector<double>(GRAV_MODEL_ORDER, 0.0)); 
std::vector<std::vector<double>> S(1, std::vector<double>(GRAV_MODEL_ORDER, 0.0)); 

Eigen::Vector3d rNedEig = Eigen::Vector3d::Zero();
Eigen::Vector3d gEcef = Eigen::Vector3d::Zero();
std::vector<double> gEci = {0.0, 0.0, 0.0};

void provideHarmonicCoefficientsJGM3(std::vector<std::vector<double>>& C, std::vector<std::vector<double>>& S) {
    C[0][0] = -0.10826360229840e-2;
    C[0][1] = -0.24140000522221e-9;
    C[0][2] =  0.15745360427672e-5;
    S[0][0] =  0.0;
    S[0][1] =  0.15430999737844e-8;
    S[0][2] = -0.90386807301869e-6;
}

void provideHarmonicCoefficientsEGM2008(std::vector<std::vector<double>>& C, std::vector<std::vector<double>>& S) {
    C.resize(GRAV_MODEL_ORDER + 1, std::vector<double>(GRAV_MODEL_ORDER + 1, 0.0));
    S.resize(GRAV_MODEL_ORDER + 1, std::vector<double>(GRAV_MODEL_ORDER + 1, 0.0));

    C[0][0] = 1.0; // Normalized coefficient for degree 0, order 0
    C[2][0] = -484.165371736e-6; // Normalized coefficient for degree 2, order 0 (J2)
    C[2][2] = 0.00000397;
    S[2][2] = 0.00000123; 
}

// Legendre Polynomial generation
double legendreP(int l, int m, double x) {
    // Compute the associated Legendre polynomial P_l^m(x)
    if (m > l) return 0.0;
    if (m == 0) {
        if (l == 0) return 1.0;
        if (l == 1) return x;
        return ((2 * l - 1) * x * legendreP(l - 1, 0, x) - (l - 1) * legendreP(l - 2, 0, x)) / l;
    } else {
        double Pnm1 = legendreP(l, m - 1, x);
        double factor = sqrt(1 - x * x);
        return factor * Pnm1;
    }
}

// Legendre Polynomial Partial Derivative
double dLegendreP(int l, int m, double x) {
    // Compute the derivative of the associated Legendre polynomial dP_l^m(x)/dx
    if (m == 0) {
        if (l == 0) return 0.0;
        if (l == 1) return 1.0;
        return (l * x * legendreP(l, 0, x) - (l + 1) * legendreP(l - 1, 0, x)) / (x * x - 1);
    } else {
        double Pnm1 = legendreP(l, m - 1, x);
        return m * x / sqrt(1 - x * x) * Pnm1;
    }
}

std::vector<double> computeGravitationalAcceleration(double r, double theta, double lambda, const std::vector<std::vector<double>>& C, const std::vector<std::vector<double>>& S) {
    double mu = G * M;
    double sin_phi = sin(M_PI / 2 - theta);
    double cos_phi = cos(M_PI / 2 - theta);

    double g_r = 0.0;
    double g_phi = 0.0;
    double g_lambda = 0.0;

    for (int l = 0; l <= GRAV_MODEL_ORDER; ++l) {
        double r_l = pow(R / r, l);

        for (int m = 0; m <= l; ++m) {
            double P_lm = legendreP(l, m, sin_phi);
            double dP_lm = dLegendreP(l, m, sin_phi);
            double cos_mlambda = cos(m * lambda);
            double sin_mlambda = sin(m * lambda);

            double C_lm = C[l][m];
            double S_lm = S[l][m];

            // Radial component
            g_r -= (mu / (r * r)) * (l + 1) * r_l * P_lm * (C_lm * cos_mlambda + S_lm * sin_mlambda);

            // Latitudinal component
            g_phi += (mu / (r * r)) * r_l * dP_lm * (C_lm * cos_mlambda + S_lm * sin_mlambda);

            // Longitudinal component
            if (m > 0) {
                g_lambda += (mu / (r * r)) * r_l * m * P_lm / cos_phi * (-C_lm * sin_mlambda + S_lm * cos_mlambda);
            }
        }
    }

    // x
    std::vector<double> g(3);
    g[0] = g_r;                     // gravitational acceleration in the r-component
    g[1] = g_phi / r;               // gravitational accelerations in latitude component
    g[2] = g_lambda / (r * cos_phi); // gravitational acceleration in longitude component

    return g;
}