#include "color-euler-method.hpp"
#include "print-color.hpp"
#include "png-writer.hpp"
#include "color-fixed-point-finite-method.hpp"
#include "shortest-path-color.hpp"
#include <cmath>
#include <iostream>
#include <string>

void random_bullet(const LCh& start, const LCh& stop) {
    double closest_theta = -M_PI;
    double closest_phi = -2.0 * M_PI;
    double closest_approach = 100.0;
    const double theta_middle = 0.5;
    const double theta_diff = 0.5;
    const double phi_middle = 1.0;
    const double phi_diff = 0.5;
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
                10000000
            );
            double min_sqr_mag = 100.0;
            for (size_t i = 0; i < out.size(); i++) {
                if (out[i].h > 4.0 * M_PI || out[i].h < -4.0 * M_PI) {
                    break;
                }
                double sqr_mag = 0.0;
                const std::array<double, 3> norm_factors{ 300.0, 800.0, 2 * M_PI };
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
        std::array<double, 3> norm_factors{ 300.0, 400.0, 2 * M_PI };
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
        actual_index = (double)(i - num_rows / 64) / (double)(62 * num_rows / 64);
        LCh interpolated_color = interpolate(out[min_index], out[0], actual_index);
        color_for_image = adobe_wide_from_XYZ(interpolated_color.to_XYZ());
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
}

void fixed_point_finite_method(const LCh& start, const LCh& stop) {
	const size_t num_rows = 2048;
	const size_t num_cols = 256;
	auto simple_shortest_path = shortest_path_approx(start, stop);
    auto out = color_finite_element_fixed_point(simple_shortest_path);
	png::PNG image_out{num_cols, num_rows};
	RGB start_color = sRGB_from_XYZ(out[0].to_XYZ());
	RGB stop_color = sRGB_from_XYZ(out.back().to_XYZ());
    for (size_t i = 0; i < num_rows / 64; i++) {
        RGB color_for_image = sRGB_from_XYZ(out[0].to_XYZ());
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
        double actual_index = (double)(i - num_rows / 64) / (double)(62 * num_rows / 64) * (double)out.size();
        size_t index = actual_index;
        XYZ lower_color = out[index].to_XYZ();
        XYZ higher_color = out[index + 1].to_XYZ();
        RGB color_for_image = sRGB_from_XYZ((lower_color + higher_color) * 0.5);
        pixel_t out_color{
	        (uint8_t)(255 * color_for_image.r),
	        (uint8_t)(255 * color_for_image.g),
	        (uint8_t)(255 * color_for_image.b),
	        255
        };
        for (size_t j = 0; j < num_cols / 4; j++) {
            image_out.at(i, j) = out_color;
        }
		lower_color = simple_shortest_path[index].to_XYZ();
		higher_color = simple_shortest_path[index + 1].to_XYZ();
		color_for_image = sRGB_from_XYZ((lower_color + higher_color) * 0.5);
		out_color = {
	        (uint8_t)(255 * color_for_image.r),
	        (uint8_t)(255 * color_for_image.g),
	        (uint8_t)(255 * color_for_image.b),
	        255
		};
        for (size_t j = num_cols / 4; j < 2 * num_cols / 4; j++) {
            image_out.at(i, j) = out_color;
        }
        actual_index = (double)(i - num_rows / 64) / (double)(62 * num_rows / 64);
        LCh interpolated_color = interpolate(out.back(), out[0], actual_index);
        color_for_image = sRGB_from_XYZ(interpolated_color.to_XYZ());
        out_color = {
	        (uint8_t)(255 * color_for_image.r),
	        (uint8_t)(255 * color_for_image.g),
	        (uint8_t)(255 * color_for_image.b),
	        255
        };
        for (size_t j = 2 * num_cols / 4; j < 3 * num_cols / 4; j++) {
            image_out.at(i, j) = out_color;
        }
		color_for_image = {
			stop_color.r * actual_index + (1.0 - actual_index) * start_color.r,
			stop_color.g * actual_index + (1.0 - actual_index) * start_color.g,
			stop_color.b * actual_index + (1.0 - actual_index) * start_color.b
		};
        out_color = {
	        (uint8_t)(255 * color_for_image.r),
	        (uint8_t)(255 * color_for_image.g),
	        (uint8_t)(255 * color_for_image.b),
	        255
        };
		for (size_t j = 3 * num_cols / 4; j < num_cols; j++) {
			image_out.at(i, j) = out_color;
		}
    }
    for (size_t i = (num_rows - num_rows / 64); i < num_rows; i++) {
        RGB color_for_image = sRGB_from_XYZ(out.back().to_XYZ());
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
    image_out.write_to_file("blue-to-orange-finite-element-hue-1.00-adjusted-distance-10000000-sRGB.png");
    full_print_color(out[0]);
	std::cout << "\n";
    full_print_color(out[out.size() / 2]);
	std::cout << "\n";
    full_print_color(out.back());
	std::cout << "\n";
}

int main() {
    LCh start{XYZ_from_sRGB({0.0, 0.0, 1.0})};
    LCh stop{XYZ_from_sRGB({1.0, 0.5, 0.0})};
    // std::swap(start, stop);
    full_print_color(start);
    std::cout << "\n";
    full_print_color({
        (start.L + stop.L) / 2.0,
        (start.C + stop.C) / 2.0,
        (start.h + stop.h) / 2.0
    });
    std::cout << "\n";
    full_print_color(stop);
    std::cout << "\n";
    fixed_point_finite_method(start, stop);
    return 0;
}
