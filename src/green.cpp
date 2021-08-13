#include "constants.hpp"
#include "XYZ.hpp"
#include "d65.h"
#include "print-color.hpp"
#include <iostream>

void green(XYZ& out, double wavelength) {
    using namespace si_constants;
    double energy = h * c / wavelength;
    if ((energy >= 1.1 * eV && energy <= 2.1 * eV)
        || (energy >= 2.5 * eV && energy <= 3.5 * eV)) {
        out = { 0., 0., 0. };
    } else {
        out = d65_rel_intensity(wavelength) * xyz_from_wavelength(wavelength);
    }
}

int main() {
    RGB out = sRGB_from_XYZ(emitted_color(d65_spectrum, green));
    std::string out_color;
    out_color.reserve(64);
    out_color += "(";
    out_color += std::to_string((int)(out.r * 255));
    out_color += ", ";
    out_color += std::to_string((int)(out.g * 255));
    out_color += ", ";
    out_color += std::to_string((int)(out.b * 255));
    out_color += ")";
    print_color(out_color.c_str(), out);
    std::cout << "\n";
    return 0;
}