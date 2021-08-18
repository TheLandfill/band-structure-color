#include "color-difference-2000.hpp"
#include "exp-by-squaring.hpp"
#include <cmath>
#include <stdexcept>

const double k_h = 1.0;
const double k_C = 1.0;
const double k_L = 1.0;

double color_distance(const LCh& c1, const LCh& c2) {
    double L1 = c1.L;
    double C1 = c1.C;
    double h1 = c1.h;
    double L2 = c2.L;
    double C2 = c2.C;
    double h2 = c2.h;
    double a1 = C1 * cos(h1);
    double b1 = C1 * sin(h1);
    double a2 = C2 * cos(h2);
    double b2 = C2 * sin(h2);
    double avg_C = C1 + (C2 - C1) / 2.0;
    double avg_C7 = avg_C;
    avg_C *= avg_C;
    avg_C7 *= avg_C;
    avg_C *= avg_C;
    avg_C7 *= avg_C;
    double G = sqrt(avg_C7 / (avg_C7 + 1808548329.0));
    double a1_prime = (1.0 + G) * a1;
    double a2_prime = (1.0 + G) * a2;
    C1 = sqrt(a1_prime * a1_prime + b1 * b1);
    C2 = sqrt(a2_prime * a2_prime + b2 * b2);
    if (a1_prime == 0.0 && b1 == 0.0) {
        h1 = 0.0;
    } else {
        h1 = atan2(b1, a1_prime);
    }
    if (a2_prime == 0.0 && b2 == 0.0) {
        h2 = 0.0;
    } else {
        h2 = atan2(b2, a2_prime);
    }
    double delta_L = L2 - L1;
    double delta_C = C2 - C1;
    double delta_h = h2 - h1;
    double mean_L = L2 + (L1 - L2) / 2.0;
    double mean_C = C2 + (C1 - C2) / 2.0;
    double mean_h = h2 + (h1 - h2) / 2.0;
    if (delta_h > M_PI) {
        delta_h -= M_PI;
    } else if (delta_h < -M_PI) {
        delta_h += M_PI;
    }
    if (C2 * C1 == 0.0) {
        delta_h = 0.0;
        mean_h *= 2.0;
    } else if (fabs(delta_h) > M_PI) {
        if (mean_h < M_PI) {
            mean_h += M_PI;
        } else {
            mean_h -= M_PI;
        }
    }
    double delta_H = 2.0 * sqrt(C1 * C2) * sin(delta_h / 2.0);
    double T = (
        -0.17 * sin(mean_h + M_PI / 3.0)
        - 0.2 * sin(4.0 * mean_h + 3.0 * M_PI / 20.0)
        + 0.24 * cos(2.0 * mean_h)
        + 0.32 * cos(3.0 * mean_h + M_PI / 30.0)
        + 1.0);

    double delta_theta_pow = -1296.0 / (25.0 * M_PI * M_PI);
    double delta_theta = 0.0;
    for (int k = -1; k <= 1; k++) {
        double offset = 2.0 * M_PI * k;
        double dist = mean_h - 55.0 * M_PI / 36.0 + offset;
        delta_theta += exp(delta_theta_pow * dist * dist);
    }
    delta_theta *= M_PI / 6.0;
    double C_to_the_7 = exp_by_squaring(mean_C, 7);
    double twenty_five_to_7 = 1808548329;
    double C7 = C_to_the_7 + twenty_five_to_7;

    double R_C = 2.0 * sqrt(C_to_the_7 / C7);

    double R_T = -sin(2.0 * delta_theta) * R_C;
    double S_C = 1.0 + 0.045 * mean_C;
    double S_H = 1.0 + 0.015 * mean_C * T;

    double L_50 = mean_L - 50.0;
    double L_50_sqr = L_50 * L_50;
    double L_50_20 = 20.0 + L_50_sqr;
    double L_50_20_sqrt = sqrt(L_50_20);
    double S_L = 1.0 + 0.015 * L_50_sqr / L_50_20_sqrt;

    double out = sqrt(
        (delta_L * delta_L / (k_L * k_L * S_L * S_L))
        + (delta_C * delta_C / (k_C * k_C * S_C * S_C))
        + (delta_H * delta_H / (k_h * k_h * S_H * S_H))
        + R_T * (delta_C * delta_H) / (k_C * S_C * k_h * S_H)
    );

    if (std::isnan(out)) {
        throw std::runtime_error("Oops. Open this in a debugger.");
    }

    return out;
}
