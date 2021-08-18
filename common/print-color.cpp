#include "print-color.hpp"
#include <iostream>

void print_color(const char * text, const RGB& color) {
	int r = color.r * 255;
	int g = color.g * 255;
	int b = color.b * 255;
	std::cout << "\x1b[38;2;" << r << ";" << g << ";" << b << "m";
	std::cout << text;
	std::cout << "\x1b[0m";
}

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
    std::cout << std::flush;
}
