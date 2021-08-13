#include "color-euler-method.hpp"
#include "print-color.hpp"
#include <cmath>
#include <iostream>
#include <string>

bool in_sRGB_gamut(double x, double y) {
    // p = p0 + (p1 - p0) * s + (p2 - p0) * t
    x -= 0.1500;
    y -= 0.0600;
    double s = 200.0 * x / 83.0 + (-500.0) * y / 747.0;
    double t = -100.0 * x / 83.0 + 4900.0 * y / 2241.0;
    return (s <= 1.0)
        && (s >= 0.0)
        && (t <= 1.0)
        && (t >= 0.0)
        && ((s + t) <= 1.0)
        && ((s + t) >= 0.0);
}

void full_print_color(const LCh& stop) {
    XYZ xyz = stop.to_XYZ();
    double total = xyz.x + xyz.y + xyz.z;
    RGB color = sRGB_from_XYZ(xyz);
    std::string out_color_message;
    out_color_message.reserve(64);
    out_color_message += "(";
    out_color_message += std::to_string(stop.L);
    out_color_message += ", ";
    out_color_message += std::to_string(stop.C);
    out_color_message += ", ";
    out_color_message += std::to_string(stop.h);
    out_color_message += ") (";
    out_color_message += std::to_string((int)(color.r * 255));
    out_color_message += ", ";
    out_color_message += std::to_string((int)(color.g * 255));
    out_color_message += ", ";
    out_color_message += std::to_string((int)(color.b * 255));
    out_color_message += ") (";
    out_color_message += std::to_string(xyz.x);
    out_color_message += ", ";
    out_color_message += std::to_string(xyz.y);
    out_color_message += ", ";
    out_color_message += std::to_string(xyz.z);
    out_color_message += ")";
    print_color(out_color_message.c_str(), color);
    xyz *= 1.0 / total;
    if (!in_sRGB_gamut(xyz.x, xyz.y)) {
        std::cout << " (OUTSIDE sRGB GAMUT)";
    }
    std::cout << "\n";
}

int main() {
    double start_L = 32.302586667249486;
    double start_a = 79.19666178930935;
    double start_b = -107.86368104495168;
    double start_C = sqrt(start_a * start_a + start_b * start_b);
    double start_h = atan2(start_b, start_a);
    double stop_L = 66.853804382266;
    double stop_a = 43.32394349110946;
    double stop_b = 73.90977076096985;
    double stop_C = sqrt(stop_a * stop_a + stop_b * stop_b);
    double stop_h = atan2(stop_b, stop_a);
    LCh start{start_L, start_C, start_h};
    LCh stop{stop_L, stop_C, stop_h};
    full_print_color(start);
    full_print_color(stop);
    full_print_color({
        (start.L + stop.L) / 2.0,
        (start.C + stop.C) / 2.0,
        (start.h + stop.h) / 2.0
    });
    double closest_theta = -M_PI;
    double closest_phi = -2.0 * M_PI;
    double closest_approach = 100.0;
    const double theta_middle = 0.5;
    const double theta_diff = 0.5;
    const double phi_middle = 1.0;
    const double phi_diff = 1.0;
    std::vector<LCh> closest_path;
    for (double theta = theta_middle - theta_diff; theta <= theta_middle + theta_diff; theta += theta_diff / 5.0) {
        #pragma omp parallel for num_threads(11)
        for (int phi_in = 0; phi_in < 11; phi_in++) {
            double phi = phi_middle - phi_diff + (double)phi_in * phi_diff / 5.0;
            auto out = out_color(
                stop,
                sin(theta * M_PI) * cos(phi * M_PI),
                sin(theta * M_PI) * sin(phi * M_PI),
                cos(theta * M_PI),
                0.0001,
                600000
            );
            double min_sqr_mag = 100.0;
            for (size_t i = 0; i < out.size(); i++) {
                if (out[i].h > 4.0 * M_PI || out[i].h < -4.0 * M_PI) {
                    break;
                }
                double sqr_mag = 0.0;
                const std::array<double, 3> norm_factors{ 100.0, 100.0, 2 * M_PI };
                for (size_t j = 0; j < 3; j++) {
                    double diff = out[i][j] - start[j];
                    diff /= norm_factors[j];
                    diff *= diff;
                    sqr_mag += diff;
                }
                if (sqr_mag < min_sqr_mag) {
                    min_sqr_mag = sqr_mag;
                }
            }
            #pragma omp critical
            {
            if (min_sqr_mag < closest_approach) {
                closest_approach = min_sqr_mag;
                closest_theta = theta;
                closest_phi = phi;
                closest_path = out;
            }
            }
            // std::cout << "Theta: " << theta << "\nPhi: " << phi << "\nClosest Approach: " << closest_approach << "\n\n";
        }
    }
    auto& out = closest_path;
    double min_sqr_mag = 10.0;
    size_t min_index = 0;
    for (size_t i = 0; i < out.size(); i++) {
        if (out[i].h > 4.0 * M_PI || out[i].h < -4.0 * M_PI) {
            break;
        }
        double sqr_mag = 0.0;
        std::array<double, 3> norm_factors{ 100.0, 100.0, 2 * M_PI };
        for (size_t j = 0; j < 3; j++) {
            double diff = out[i][j] - start[j];
            diff /= norm_factors[j];
            diff *= diff;
            sqr_mag += diff;
        }
        if (sqr_mag < min_sqr_mag) {
            min_sqr_mag = sqr_mag;
            min_index = i;
        }
    }
    //for (size_t i = 0; i < out.size(); i += 600) {
    //    full_print_color(out[i]);
    //}
    //for (size_t i = 0; i < out.size(); i += 600) {
    //    XYZ color = out[i].to_XYZ();
    //    double norm = color.x + color.y + color.z;
    //    color *= 1.0 / norm;
    //    std::cout << color.x << "," << color.y << "\n";
    //}
    std::cout << min_index << ": ";
    full_print_color(out[min_index]);
    std::cout << min_index / 2 << ": ";
    full_print_color(out[min_index / 2]);
    std::cout << "000000: ";
    full_print_color(out[0]);
    std::cout << "Theta: " << closest_theta << "\nPhi: " << closest_phi << "\nApproach: " << closest_approach << "\n";
    std::cout << "Initial Velocity: (" <<  sin(closest_theta * M_PI) * cos(closest_phi * M_PI) << ", "
        << sin(closest_theta * M_PI) * sin(closest_phi * M_PI) << ", "
        << cos(closest_theta * M_PI) << ")\n";
    return 0;
}