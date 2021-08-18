#include "LCh.hpp"
#include "exp-by-squaring.hpp"
#include <stdexcept>
#include <cmath>

LCh::LCh() : L(0.0), C(0.0), h(0.0) {}

LCh::LCh(double luma, double chroma, double hue) : L(luma), C(chroma), h(hue) {}

double lab_from_xyz_f(double t);

LCh::LCh(XYZ in) {
    static const double d65_X_n = 95.0489 / 100.0;
    static const double d65_Y_n = 100.0 / 100.0;
    static const double d65_Z_n = 108.8840 / 100.0;
    double x_f = lab_from_xyz_f(in.x / d65_X_n);
    double y_f = lab_from_xyz_f(in.y / d65_Y_n);
    double z_f = lab_from_xyz_f(in.z / d65_Z_n);
    L = 116.0 * y_f - 16.0;
    double a = 500.0 * (x_f - y_f);
    double b = 200.0 * (y_f - z_f);
    C = sqrt(a * a + b * b);
    double G = 0.5 * (1.0 -
        sqrt(
            exp_by_squaring(C, 7) / (exp_by_squaring(C, 7) + exp_by_squaring(25.0, 7))
        )
    );
    a *= (1.0 + G);
    C = sqrt(a * a + b * b);
    h = atan2(b, a);
}

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
    double G = 0.5 * (1.0 -
        sqrt(
            exp_by_squaring(C, 7) / (exp_by_squaring(C, 7) + exp_by_squaring(25.0, 7))
        )
    );
    a /= (1.0 + G);
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

LCh close_interpolate(const LCh& a, const LCh& b, double t) {
	return {
        a.L * t + (1.0 - t) * b.L,
        a.C * t + (1.0 - t) * b.C,
		a.h * t + (1.0 - t) * b.h
	};
}

LCh interpolate(const LCh& a, const LCh& b, double t) {
    double h1 = a.h;
    double h2 = b.h;
    if (fabs(h1 - h2) > M_PI) {
        if (h1 > h2) {
            h1 -= 2.0 * M_PI;
        } else {
            h2 -= 2.0 * M_PI;
        }
    }
    double hue = h1 * t + (1.0 - t) * h2;
    if (hue >= M_PI * 2.0) {
        hue -= 2.0 * M_PI;
    } else if (hue <= 0.0) {
        hue += 2.0 * M_PI;
    }
    return {
        a.L * t + (1.0 - t) * b.L,
        a.C * t + (1.0 - t) * b.C,
        hue
    };
}
