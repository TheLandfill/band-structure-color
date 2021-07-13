#include "blackbody.hpp"
#include "constants.hpp"
#include "XYZ.hpp"
#include <cmath>

Black_Body::Black_Body(double T) : temperature(T) {}

// Units are W / m^2 / m / steradian
double Black_Body::operator()(double wavelength) {
    using namespace si_constants;
    double numerator = 2.0 * h * c * c;
    double wavelength_pow = wavelength * wavelength;
    wavelength_pow *= wavelength_pow;
    wavelength_pow *= wavelength;
    double exponent = h * c;
    exponent /= wavelength;
    exponent /= kb;
    exponent /= temperature;
    double boltzman_denominator = exp(exponent) - 1.0;
    double result = numerator / wavelength_pow;
    result /= boltzman_denominator;
    return result;
}

void Black_Body::operator()(double& out, double wavelength) {
    out = (*this)(wavelength);
}

void Black_Body::operator()(XYZ& out, double wavelength) {
    out = (*this)(wavelength) * xyz_from_wavelength(wavelength);
}