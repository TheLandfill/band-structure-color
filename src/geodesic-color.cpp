#include "color-euler-method.hpp"
#include "print-color.hpp"
#include "png-writer.hpp"
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
    RGB color = adobe_wide_from_XYZ(xyz);
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
    LCh start{XYZ_from_adobe_wide({0.918, 0.451, 0.195})};
    // LCh stop{20.346538, 130.120806, -0.755002};
    LCh stop{XYZ_from_adobe_wide({0.10, 0.10, 1.0})};
    // LCh swap_color = start;
    // start = stop;
    // stop = swap_color;
    // stop.L = start.L;
    full_print_color(start);
    full_print_color({
        (start.L + stop.L) / 2.0,
        (start.C + stop.C) / 2.0,
        (start.h + stop.h) / 2.0
    });
    full_print_color(stop);
    double closest_theta = -M_PI;
    double closest_phi = -2.0 * M_PI;
    double closest_approach = 100.0;
    const double theta_middle = 0.508125;
    const double theta_diff = 0.5 / 20.0 / 20.0 / 20.0;
    const double phi_middle = 1.0475;
    const double phi_diff = 1.0 / 20.0 / 20.0 / 20.0;
    std::vector<LCh> closest_path;
    for (double theta = theta_middle - theta_diff; theta <= theta_middle + theta_diff; theta += theta_diff / 20.0) {
        #pragma omp parallel for num_threads(11)
        for (int phi_in = 0; phi_in < 41; phi_in++) {
            double phi = phi_middle - phi_diff + (double)phi_in * phi_diff / 20.0;
            auto out = out_color(
                start,
                sin(theta * M_PI) * cos(phi * M_PI),
                sin(theta * M_PI) * sin(phi * M_PI),
                cos(theta * M_PI),
                0.0001,
                1000000
            );
            double min_sqr_mag = 100.0;
            for (size_t i = 0; i < out.size(); i++) {
                if (out[i].h > 4.0 * M_PI || out[i].h < -4.0 * M_PI) {
                    break;
                }
                double sqr_mag = 0.0;
                const std::array<double, 3> norm_factors{ 300.0, 290.0, 2 * M_PI };
                for (size_t j = 0; j < 3; j++) {
                    double diff = out[i][j] - stop[j];
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
        std::array<double, 3> norm_factors{ 300.0, 290.0, 2 * M_PI };
        for (size_t j = 0; j < 3; j++) {
            double diff = out[i][j] - stop[j];
            diff /= norm_factors[j];
            diff *= diff;
            sqr_mag += diff;
        }
        if (sqr_mag < min_sqr_mag) {
            min_sqr_mag = sqr_mag;
            min_index = i;
        }
    }
    // for (size_t i = 0; i < out.size(); i += 4000) {
    //     full_print_color(out[i]);
    // }
    // for (size_t i = 0; i < out.size(); i += 600) {
    //     XYZ color = out[i].to_XYZ();
    //     double norm = color.x + color.y + color.z;
    //     color *= 1.0 / norm;
    //     std::cout << color.x << "," << color.y << "\n";
    // }
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
	const size_t num_rows = 2048;
	const size_t num_cols = 256;
	png::PNG image_out{num_cols, num_rows};
    for (size_t i = 0; i < num_rows / 64; i++) {
        RGB color_for_image = adobe_wide_from_XYZ(out[0].to_XYZ());
        pixel_t out_color{
	        (uint8_t)(255 * color_for_image.r), 
	        (uint8_t)(255 * color_for_image.g),
	        (uint8_t)(255 * color_for_image.b),
	        255
        };
        for (size_t j = 0; j < num_cols; j++) {
            image_out.at(i, j) = out_color;
        }
    }
    for (size_t i = num_rows / 64; i < (num_rows - num_rows / 64); i++) {
        double actual_index = (double)(i - num_rows / 64) / (double)(62 * num_rows / 64) * (double)min_index;
        size_t index = actual_index;
        XYZ lower_color = out[index].to_XYZ();
        XYZ higher_color = out[index + 1].to_XYZ();
        RGB color_for_image = adobe_wide_from_XYZ((lower_color + higher_color) * 0.5);
        pixel_t out_color{
	        (uint8_t)(255 * color_for_image.r), 
	        (uint8_t)(255 * color_for_image.g),
	        (uint8_t)(255 * color_for_image.b),
	        255
        };
        for (size_t j = 0; j < num_cols; j++) {
            image_out.at(i, j) = out_color;
        }
        actual_index = (double)i / (double)num_rows;
        lower_color = out[0].to_XYZ() * (1.0 - actual_index);
        higher_color = out[min_index].to_XYZ() * (actual_index);
        color_for_image = adobe_wide_from_XYZ((lower_color + higher_color));
        out_color = {
	        (uint8_t)(255 * color_for_image.r), 
	        (uint8_t)(255 * color_for_image.g),
	        (uint8_t)(255 * color_for_image.b),
	        255
        };
        for (size_t j = num_cols / 2; j < num_cols; j++) {
            image_out.at(i, j) = out_color;
        }
    }
    for (size_t i = (num_rows - num_rows / 64); i < num_rows; i++) {
        RGB color_for_image = adobe_wide_from_XYZ(out[min_index].to_XYZ());
        pixel_t out_color{
	        (uint8_t)(255 * color_for_image.r), 
	        (uint8_t)(255 * color_for_image.g),
	        (uint8_t)(255 * color_for_image.b),
	        255
        };
        for (size_t j = 0; j < num_cols; j++) {
            image_out.at(i, j) = out_color;
        }
    }
    image_out.write_to_file("blue-to-orange-fixed-maybe.png");
    return 0;
}