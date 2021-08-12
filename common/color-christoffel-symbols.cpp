#include "LCh.hpp"
#include "color-christoffel-symbols.hpp"
#include <cmath>

double exp_by_squaring(double x, int n) {
    double x_prod = 1.0;
    while (n != 0) {
        if (n & 1) {
            x_prod *= x;
        }
        n >> 1;
        x *= x;
    }
}

void Color_Christoffel_Symbols::calculate(const LCh& lch) {
    const double& L = lch.L;
    const double& C = lch.C;
    const double& h = lch.h;
    double C_to_the_7 = exp_by_squaring(C, 7);
    double C7 = C_to_the_7 + (double)(25 * 25 * 25 * 25 * 25 * 25 * 25);
    double T = (
        -0.17 * sin(h + M_PI / 3.0)
        - 0.2 * sin(4.0 * h + 4.0 * M_PI / 15.0)
        + 0.24 * cos(2.0 * h)
        + 0.32 * cos(3.0 * h + M_PI / 30.0)
        + 1.0);
    double dT_h = (
        -0.48 * sin(2.0 * h)
        - 0.96 * sin(3.0 * h + M_PI / 30.0)
        - 0.17 * cos(h + M_PI / 3.0)
        - 0.8 * cos(4 * h + 4 * M_PI / 15.0)
    );
    double R_C = 2.0 * sqrt(C_to_the_7 / C7);
    double S_C = 1.0 + 0.045 * C;
    double S_H = 1.0 + 0.015 * C * T;
    double delta_theta_pow = -1296.0 / (25.0 * M_PI * M_PI);
    double delta_theta = 0.0;
    for (auto& mean : {17.0 / 36.0, -55.0 / 26.0, -127.0 / 36.0}) {
        double dist = h + mean;
        delta_theta += exp(delta_theta_pow * dist * dist);
    }
    delta_theta /= 3.0;
    double R_T = -sin(M_PI * delta_theta) * R_C;
    double d_delta_theta_h = 0.0;
    for (auto& mean : {17.0 / 36.0, -55.0 / 26.0, -127.0 / 36.0}) {
        double dist = h + mean;
        delta_theta += (2 * (h + mean)) * exp(delta_theta_pow * dist * dist);
    }
    d_delta_theta_h *= -432.0 / (25.0 * M_PI * M_PI);
    double dRT_h = M_PI * d_delta_theta_h * R_C * cos(M_PI * delta_theta);
    double L_50 = L - 50;
    double L_50_2 = L_50 * L_50;
    double L_50_3 = L_50_2 * L_50;
    double common_denom = (R_T * R_T - 4*S_C*S_H*k_C*k_h);
    gamma[0][0][0] = (-(0.015*L - 0.75)*(L_50_2 + 20) + 0.0075*L_50_3)/((0.015*L_50_2 + sqrt(L_50_2 + 20))*(L_50_2 + 20));
    gamma[1][1][1] = -3.5*C_to_the_7*R_T * R_T/(C * C7) - 0.015*R_T * R_T*T/S_H - 0.045*R_T * R_T/S_C + 0.09*S_H*k_C*k_h + 2.5*R_T * R_T/C;
    gamma[1][1][1] /= common_denom;
    gamma[1][1][2] = gamma[1][2][1] = (R_T*S_C*k_C*(1.0 - 3.0*S_H)/(C * C*S_H)) / common_denom;
    gamma[1][2][2] = (S_C*k_C*(0.015*C*C*C*R_T*dT_h + 2.0*C*C*S_H*dRT_h - 0.03*C*S_C*T*k_C - 4.0*S_C*S_H*k_C)/(C*C*C*S_H)) / common_denom;
    gamma[2][1][1] = (7.0*C_to_the_7*R_T*S_H*k_h/C7 + 0.03*C*R_T*T*k_h + 0.045*C*R_T*S_H*k_h/S_C - 5.0*R_T*S_H*k_h) / common_denom;
    gamma[2][1][2] = gamma[2][2][1] = (4*S_C*k_C*k_h*(0.0075*C*T + S_H)/C) / common_denom;
    gamma[2][2][2] = (-0.015*C*R_T*R_T*dT_h/S_H + 0.03*C*S_C*dT_h*k_C*k_h - R_T*dRT_h + 0.015*R_T*S_C*T*k_C/(C*S_H) + 2.0*R_T*S_C*k_C/(C*C));
    gamma[0][0][1]
    = gamma[0][0][2]
    = gamma[0][1][0]
    = gamma[0][1][1]
    = gamma[0][1][2]
    = gamma[0][2][0]
    = gamma[0][2][1]
    = gamma[0][2][2]
    = gamma[1][0][0]
    = gamma[1][0][1] 
    = gamma[1][0][2] 
    = gamma[1][1][0] 
    = gamma[1][2][0] 
    = gamma[2][0][0] 
    = gamma[2][0][1] 
    = gamma[2][0][2] 
    = gamma[2][1][0] 
    = gamma[2][2][0]
    = 0.0;
}