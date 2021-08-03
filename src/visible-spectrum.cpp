#include "XYZ.hpp"
#include "template-integration.hpp"
#include "blackbody.hpp"
#include "d65.h"
#include "constants.hpp"
#include "png-writer.hpp"
#include <string>
#include <iostream>

class Small_Color_Range {
public:
	Small_Color_Range(double energy) :
        min_energy((energy - 0.01) * si_constants::eV),
        max_energy((energy + 0.01) * si_constants::eV)
        {}
	void operator()(XYZ& out, double wavelength) {
		using namespace si_constants;
        double energy = h * c / wavelength;
		if (energy < min_energy || energy > max_energy) {
			out = { 0., 0., 0. };
		} else {
			out = d65_rel_intensity(wavelength) * xyz_from_wavelength(wavelength);
		}
	}
private:
	double min_energy;
	double max_energy;
};


int main() {
	const size_t num_rows = 2048;
	const size_t num_cols = 256;
	png::PNG out{num_cols, num_rows};
	const double start = 1.75;
	const double end = 3.25;
	const double num_steps = num_rows;
	const double step_size = (end - start) / num_steps;
	size_t cur_row = 0;
	for (double i = start; i <= end; i += step_size) {
		Small_Color_Range sc{i};
		XYZ output = emitted_color(d65_spectrum, sc);
		RGB s_out = sRGB_from_XYZ(output);
        pixel_t out_color{0, 0, 0, 255};
        out_color = {
			(uint8_t)(255 * s_out.r), 
			(uint8_t)(255 * s_out.g),
			(uint8_t)(255 * s_out.b),
			255
        };
		for (size_t cur_col = 0; cur_col < num_cols; cur_col++) {
			out.at(cur_row, cur_col) = out_color;
		}
		cur_row++;
	}
	out.write_to_file("./visible-spectrum.png");
	return 0;
}