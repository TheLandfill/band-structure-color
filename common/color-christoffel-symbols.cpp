#include "LCh.hpp"
#include "color-christoffel-symbols.hpp"
#include <cmath>

double exp_by_squaring(double x, int n) {
    double x_prod = 1.0;
    while (n != 0) {
        if (n & 1) {
            x_prod *= x;
        }
        n >>= 1;
        x *= x;
    }
    return x_prod;
}

void Color_Christoffel_Symbols::calculate(const LCh& lch) {
    const double& L = lch.L;
    const double& C = lch.C;
    const double& h = lch.h;
    // Calculating intermediate variables for a bit
    double T = (
        -0.17 * sin(h + M_PI / 3.0)
        - 0.2 * sin(4.0 * h + 3.0 * M_PI / 20.0)
        + 0.24 * cos(2.0 * h)
        + 0.32 * cos(3.0 * h + M_PI / 30.0)
        + 1.0);
    double dT_h = (
        -0.48 * sin(2.0 * h)
        - 0.96 * sin(3.0 * h + M_PI / 30.0)
        - 0.17 * cos(h + M_PI / 3.0)
        - 0.8 * cos(4 * h + 3.0 * M_PI / 20.0)
    );

    double C_to_the_7 = exp_by_squaring(C, 7);
    double twenty_five_to_7 = 1808548329;
    double C7 = C_to_the_7 + twenty_five_to_7;

    double R_C = 2.0 * sqrt(C_to_the_7 / C7);
    double dR_C = R_C * (7.0 / 2.0) / (C * (C_to_the_7 / twenty_five_to_7 + 1.0));

    double S_C = 1.0 + 0.045 * C;
    double dS_C = 0.045;

    double S_H = 1.0 + 0.015 * C * T;
    double dS_H_C = 0.015 * T;
    double dS_H_h = 0.015 * C * dT_h;

    double delta_theta_pow = -1296.0 / (25.0 * M_PI * M_PI);

    double delta_theta = 0.0;
    for (int k = -1; k <= 1; k++) {
        double offset = 2.0 * M_PI * k;
        double dist = h - 55.0 * M_PI / 36.0 + offset;
        delta_theta += exp(delta_theta_pow * dist * dist);
    }
    delta_theta *= M_PI / 6.0;

    double R_T = -sin(2.0 * delta_theta) * R_C;

    double dRT_C = -sin(2.0 * delta_theta) * dR_C;

    double d_delta_theta_h = 0.0;
    for (int k = -1; k <= 1; k++) {
        double offset = 2.0 * M_PI * (double)k;
        double dist = h - 55.0 * M_PI / 36.0 + offset;
        d_delta_theta_h += exp(delta_theta_pow * dist * dist) * 2.0 * dist * delta_theta_pow;
    }
    d_delta_theta_h *= M_PI / 6.0;

    double dRT_h = -2.0 * d_delta_theta_h * R_C * cos(2.0 * delta_theta);

    double L_50 = L - 50.0;
    double L_50_sqr = L_50 * L_50;
    double L_50_20 = 20.0 + L_50_sqr;
    double L_50_20_sqrt = sqrt(L_50_20);
    double S_L = 1.0 + 0.015 * L_50_sqr / L_50_20_sqrt;
    double dS_L = (3.0 / 200.0) * L_50 * (L_50_20 + 20.0) / (L_50_20 * L_50_20_sqrt);

    double common = R_T * R_T - 4.0;

    // Calculating the Christoffel symbols

    // 1
    gamma[0][0][0] = dS_L / S_L;
    // 2, 3
    gamma[2][1][2] = gamma[2][2][1] = 4.0 * (C * dS_H_C - S_H) / (C * common * S_H);
    // 4, 5
    gamma[1][1][2] = gamma[1][2][1] = -(C * k_C * R_T * S_C) / (2.0 * k_h * S_H) * gamma[2][1][2];
    // 6
    gamma[1][1][1] = (
        C * dRT_C * R_T * S_C * S_H
        - common * C * dS_C * S_H
        - C * R_T * R_T * S_C * dS_H_C
        + R_T * R_T * S_C * S_H
    ) / (
        C * common * S_C * S_H
    );
    // 7
    gamma[1][2][2] = C * k_C * S_C * (
        4.0 * k_C * S_C * S_H
        - 2.0 * dRT_h * k_h * S_H * S_H
        - 4.0 * C * dS_H_C * k_C * S_C
    ) / (
        k_h * k_h * common * S_H * S_H * S_H
    );
    // 8
    gamma[2][1][1] = 2.0 * k_h * (
        C * R_T * dS_H_C
        - C * dRT_h * S_H
        - R_T * S_H
    ) / (
        C * C * k_C * common * S_C
    );
    // 9
    gamma[2][2][2] = (
        2.0 * C * dS_H_C * k_C * R_T * S_C
        + dRT_h * k_h * R_T * S_H * S_H
        - dS_H_h * k_h * R_T * R_T * S_H
        + 4.0 * dS_H_h * k_h * S_H
        - 2.0 * k_C * R_T * S_C * S_H
    ) / (
        k_h * common * S_H * S_H
    ); 

    // All of these are zero

    // 8 for [0][a][b] where a b != 0
    // 17
    gamma[0][0][1]
    = gamma[0][0][2]
    = gamma[0][1][0]
    = gamma[0][1][1]
    = gamma[0][1][2]
    = gamma[0][2][0]
    = gamma[0][2][1]
    = gamma[0][2][2]
    // 5 for [1][a][b] where a b == 0
    // 22
    = gamma[1][0][0]
    = gamma[1][0][1] 
    = gamma[1][0][2] 
    = gamma[1][1][0] 
    = gamma[1][2][0] 
    // 5 for [1][a][b] where a b == 0
    // 27
    = gamma[2][0][0] 
    = gamma[2][0][1] 
    = gamma[2][0][2] 
    = gamma[2][1][0] 
    = gamma[2][2][0]
    = 0.0;
}