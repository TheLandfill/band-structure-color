#include "LCh.hpp"
#include <stdexcept>
#include <cmath>

LCh::LCh(double luma, double chroma, double hue) : L(luma), C(chroma), h(hue) {}

double lab_from_xyz_f(double t) {
    double delta = 6.0 / 29.0;
    if (t > delta * delta * delta) {
        return pow(t, 1.0 / 3.0);
    } else {
        return t / (3.0 * delta * delta) + 4.0 / 29.0;
    }
}

double xyz_from_lab_f(double t) {
    double delta = 6.0 / 29.0;
    if (t > delta) {
        return t * t * t;
    } else {
        return 3.0 * delta * delta * (t - 4.0 / 29.0);
    }
}

XYZ LCh::to_XYZ() const {
    double a = C * cos(h);
    double b = C * sin(h);
    static const double d65_X_n = 95.0489 / 100.0;
    static const double d65_Y_n = 100.0 / 100.0;
    static const double d65_Z_n = 108.8840 / 100.0;
    double term_L = (L + 16.0) / 116.0;
    return {
        d65_X_n * xyz_from_lab_f(term_L + a / 500.0),
        d65_Y_n * xyz_from_lab_f(term_L),
        d65_Z_n * xyz_from_lab_f(term_L - b / 200.0)
    };
}

double& LCh::operator[](size_t i) {
    switch(i) {
    case 0:
        return L;
    case 1:
        return C;
    case 2:
        return h;
    default:
        throw std::out_of_range("ERROR in LCh index\n");
    }
}

double LCh::operator[](size_t i) const {
    switch(i) {
    case 0:
        return L;
    case 1:
        return C;
    case 2:
        return h;
    default:
        throw std::out_of_range("ERROR in LCh index\n");
    }
}