#pragma once

#ifdef __cplusplus
extern "C" {
#endif
struct XYZ {
	double x, y, z;
};

struct RGB {
	double r, g, b;
};

struct RGB XYZ_to_RGB(const struct XYZ * in);

struct XYZ xyz_from_wavelength(double lambda);

#ifdef __cplusplus
}
#endif