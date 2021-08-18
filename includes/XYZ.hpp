#pragma once
#include "XYZ.h"
#include <functional>

XYZ& operator+=(XYZ& a, const XYZ& b);
XYZ& operator*=(XYZ& a, const double scale);
XYZ operator*(const double scale, const XYZ& a);
XYZ operator*(const XYZ& a, const double scale);
XYZ operator+(const XYZ& a, const XYZ& b);

typedef void (*N_func)(double&, double);
typedef std::function<void(XYZ&, double)> XYZ_func;

XYZ emitted_color(N_func, XYZ_func);
void d65_spectrum(double& out, double wavelength);
void d65_spectrum(XYZ& out, double wavelength);

RGB sRGB_from_XYZ(const XYZ& in);
XYZ XYZ_from_sRGB(const RGB& in);
RGB adobe_wide_from_XYZ(const XYZ& in);
XYZ XYZ_from_adobe_wide(RGB in);
