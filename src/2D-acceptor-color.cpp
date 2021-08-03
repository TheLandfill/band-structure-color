#include "XYZ.hpp"
#include "template-integration.hpp"
#include "blackbody.hpp"
#include "d65.h"
#include "constants.hpp"
#include "png-writer.hpp"
#include <string>
#include <iostream>

// Top of band is zero
class Bottom_and_Acceptor_Level_Gap_Function {
public:
	Bottom_and_Acceptor_Level_Gap_Function(double bottom_of_band, double acceptor_level) :
		min_energy(acceptor_level * si_constants::eV),
		max_energy((acceptor_level - bottom_of_band) * si_constants::eV)
	{}
	void operator()(XYZ& out, double wavelength) {
		using namespace si_constants;
		double energy = h * c / wavelength;
		if ((energy >= min_energy) && (energy <= max_energy)) {
			out = { 0., 0., 0. };
		} else {
			out = d65_rel_intensity(wavelength) * xyz_from_wavelength(wavelength);
		}
	}
	void print_energy_range() {
		std::cout << "Min: " << min_energy << "\nMax: " << max_energy << "\n";
	}
private:
	double min_energy;
	double max_energy;
};

int main() {
	const size_t num_rows = 2048;
	const size_t num_cols = 2048;
	png::PNG out{num_cols, num_rows};
	const double bottom_start = 0.0;
	const double bottom_end = -3.0;
	const double bottom_step = (bottom_end - bottom_start) / num_rows;
	const double acceptor_start = 0.0;
	const double acceptor_end = 2.75;
	const double acceptor_step = (acceptor_end - acceptor_start) / num_cols;
	#pragma omp parallel for num_threads(12)
	for (size_t i = 0; i < num_rows; i++) {
		double bottom = bottom_start + bottom_step * i;
		double acceptor = acceptor_start;
		for (size_t j = 0; j < num_cols; j++) {
			Bottom_and_Acceptor_Level_Gap_Function f{bottom, acceptor};
			XYZ output = emitted_color(d65_spectrum, f);
			RGB s_out = sRGB_from_XYZ(output);
			out.at(i, j) = {
				(uint8_t)(255 * s_out.r), 
				(uint8_t)(255 * s_out.g),
				(uint8_t)(255 * s_out.b),
				255
			};
			acceptor += acceptor_step;
		}
	}
	out.write_to_file("./acceptor-energy-levels-updated.png");
	return 0;
}
