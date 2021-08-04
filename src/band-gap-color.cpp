#include "XYZ.hpp"
#include "template-integration.hpp"
#include "blackbody.hpp"
#include "d65.h"
#include "constants.hpp"
#include "png-writer.hpp"
#include <string>
#include <iostream>

class Simple_Gap_Function {
public:
	Simple_Gap_Function(double energy) : min_energy(energy) {}
	void operator()(XYZ& out, double wavelength) {
		using namespace si_constants;
		if ((h * c / wavelength) <= min_energy * eV) {
			out = { 0., 0., 0. };
		} else {
			out = d65_rel_intensity(wavelength) * xyz_from_wavelength(wavelength);
		}
	}
private:
	double min_energy;
};

int main() {
	const size_t num_rows = 2048;
	const size_t num_cols = 256;
	png::PNG out{num_cols, num_rows};
	std::string out_text;
	out_text.reserve(128);
	const double start = 1.75;
	const double end = 3.25;
	const double num_steps = num_rows;
	const double step_size = (end - start) / num_steps;
	size_t cur_row = 0;
	for (double i = start; i <= end; i += step_size) {
		out_text.clear();
		Simple_Gap_Function sg{i};
		XYZ output = emitted_color(d65_spectrum, sg);
		RGB s_out = sRGB_from_XYZ(output);
		for (size_t cur_col = 0; cur_col < num_cols; cur_col++) {
			out.at(cur_row, cur_col) = {
				(uint8_t)(255 * s_out.r), 
				(uint8_t)(255 * s_out.g),
				(uint8_t)(255 * s_out.b),
				255
			};
		}
		cur_row++;
	}
	out.write_to_file("./reversed-band-gap.png");
	return 0;
}
