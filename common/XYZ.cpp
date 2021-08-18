#include "XYZ.hpp"
#include "template-integration.hpp"
#include "d65.h"
#include <cmath>
#include <algorithm>
#include <iostream>

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
    //output *= 1.0 / (output.x + output.y + output.z);
    // std::cout << output.x << ", " << output.y << ", " << output.z << ":\t";
    //output *= 1.0 / output.y;
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

void reverse_gamma_correction(double& d) {
	if (d <= 0.4045) {
		d /= 12.92;
	} else {
		static const double a = 0.055;
		d = pow((d + a) / (1.0 + a), 2.4);
	}
}

XYZ XYZ_from_sRGB(const RGB& in) {
	RGB copy = in;
	reverse_gamma_correction(copy.r);
	reverse_gamma_correction(copy.g);
	reverse_gamma_correction(copy.b);
	return {
		0.4124 * copy.r + 0.3576 * copy.g + 0.1805 * copy.b,
		0.2126 * copy.r + 0.7152 * copy.g + 0.0722 * copy.b,
		0.0193 * copy.r + 0.1192 * copy.g + 0.9505 * copy.b
	};
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

RGB adobe_wide_from_XYZ(const XYZ& in) {
    RGB out {
         1.4628067 * in.x - 0.1840623 * in.y - 0.2743606 * in.z,
        -0.5217933 * in.x + 1.4472381 * in.y + 0.0677227 * in.z,
         0.0349342 * in.x - 0.0968930 * in.y + 1.2884099 * in.z
    };
    out.r = std::clamp(out.r, 0.0, 1.0);
    out.g = std::clamp(out.g, 0.0, 1.0);
    out.b = std::clamp(out.b, 0.0, 1.0);
    out.r = pow(out.r, 1.0 / 2.19921875);
    out.g = pow(out.g, 1.0 / 2.19921875);
    out.b = pow(out.b, 1.0 / 2.19921875);
    out.r = std::clamp(out.r, 0.0, 1.0);
    out.g = std::clamp(out.g, 0.0, 1.0);
    out.b = std::clamp(out.b, 0.0, 1.0);
    return out;
}

XYZ XYZ_from_adobe_wide(RGB in) {
    in.r = pow(in.r, 2.19921875);
    in.g = pow(in.g, 2.19921875);
    in.b = pow(in.b, 2.19921875);
    XYZ out = {
        0.7161046 * in.r + 0.1009296 * in.g + 0.1471858 * in.b,
        0.2581874 * in.r + 0.7249378 * in.g + 0.0168748 * in.b,
        0.0000000 * in.r + 0.0517813 * in.g + 0.7734287 * in.b
    };
    return out;
}
