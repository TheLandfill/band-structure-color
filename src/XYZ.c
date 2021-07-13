#include "XYZ.h"
#include <math.h>

struct RGB XYZ_to_RGB(const struct XYZ* in) {
	struct RGB out;
	out.r =  0.41847   * in->x - 0.15866   * in->y - 0.082835 * in->z;
	out.g = -0.091169  * in->x + 0.25243   * in->y + 0.015708 * in->z;
	out.b =  0.0009209 * in->x - 0.0025498 * in->y + 0.17860  * in->z;
	return out;
};

double gaussian(double x, double scale, double mean, double std_left, double std_right) {
	double t = (x - mean) / (x < mean ? std_left : std_right);
	return scale * exp(-(t * t) / 2.0);
}

struct XYZ xyz_from_wavelength(double lambda) {
	struct XYZ color;
	color.x =	gaussian(lambda,  1.056, 599.8e-9, 37.9e-9, 31.0e-9) +
				gaussian(lambda,  0.362, 442.0e-9, 16.0e-9, 26.7e-9) +
				gaussian(lambda, -0.065, 501.1e-9, 20.4e-9, 26.2e-9);
	color.y =	gaussian(lambda,  0.821, 568.8e-9, 46.9e-9, 40.5e-9) +
				gaussian(lambda,  0.286, 530.9e-9, 16.3e-9, 31.1e-9);
	color.z =	gaussian(lambda,  1.217, 437.0e-9, 11.8e-9, 36.0e-9) +
				gaussian(lambda,  0.681, 459.0e-9, 26.0e-9, 13.8e-9);
	return color;
}
