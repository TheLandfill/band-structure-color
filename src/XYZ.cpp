#include "XYZ.hpp"
#include "template-integration.hpp"
#include "d65.h"
#include <cmath>
#include <algorithm>

XYZ& operator+=(XYZ& a, const XYZ& b) {
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
	return a;
}

XYZ& operator*=(XYZ& a, const double scale) {
	a.x *= scale;
	a.y *= scale;
	a.z *= scale;
	return a;
}

XYZ operator*(const double scale, const XYZ& a) {
	return { a.x * scale, a.y * scale, a.z * scale };
}

XYZ operator*(const XYZ& a, const double scale) {
	return { a.x * scale, a.y * scale, a.z * scale };
}

XYZ operator+(const XYZ& a, const XYZ& b) {
	return { a.x + b.x, a.y + b.y, a.z + b.z };
}

XYZ emitted_color(N_func N, XYZ_func xyz) {
    Simpson_Points spectrum{300e-9, 830e-9, 2048};
    Simpson_Sums<double> spectrum_sum;
    Simpson_Sums<XYZ> full_spectrum_sum;
    double n = simpson_rule(N, spectrum, spectrum_sum);
    XYZ output = simpson_rule(xyz, spectrum, full_spectrum_sum);
    output *= 1.0 / n;
    output *= 1.0 / (output.x + output.y + output.z);
    output *= 1.0 / output.y;
    return output;
}

void d65_spectrum(double& out, double wavelength) {
    out = d65_rel_intensity(wavelength) * xyz_from_wavelength(wavelength).y;
}

void d65_spectrum(XYZ& out, double wavelength) {
    out = d65_rel_intensity(wavelength) * xyz_from_wavelength(wavelength);
}

void gamma_correction(double& d) {
    if (d <= 0.0031308) {
        d *= 12.92;
    } else {
        d = pow(d, 1.0/2.4);
        d *= 1.055;
        d -= 0.055;
    }
}

RGB sRGB_from_XYZ(const XYZ& in) {
    RGB out{
         3.2406 * in.x - 1.5372 * in.y - 0.4986 * in.z,
        -0.9689 * in.x + 1.8758 * in.y + 0.0415 * in.z,
         0.0557 * in.x - 0.2040 * in.y + 1.0570 * in.z
    };
    gamma_correction(out.r);
    gamma_correction(out.g);
    gamma_correction(out.b);
    out.r = std::clamp(out.r, 0.0, 1.0);
    out.g = std::clamp(out.g, 0.0, 1.0);
    out.b = std::clamp(out.b, 0.0, 1.0);
    return out;
}