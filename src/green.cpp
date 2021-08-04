#include "constants.hpp"
#include "XYZ.hpp"

class Small_Color_Range {
public:
	Small_Color_Range(double energy) :
        min_energy((energy - 7.0 / 32.0) * si_constants::eV),
        max_energy((energy + 7.0 / 32.0) * si_constants::eV)
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
    
    return 0;
}