#include "XYZ.hpp"
#include "template-integration.hpp"
#include "blackbody.hpp"
#include "d65.h"
#include "constants.hpp"
#include <string>
#include <iostream>

void sample_func(double& value, double x) {
	value = 1;
}

void sunlight_func(double& a, double b) {
	static Black_Body sunlight{6504.0};
	return sunlight(a, b);
}

void sunlight_func(XYZ& a, double b) {
	static Black_Body sunlight{6504.0};
	return sunlight(a, b);
}

void rust_func(XYZ& out, double wavelength) {
	using namespace si_constants;
	if ((h * c / wavelength) >= 2.1 * eV) {
		out = { 0., 0., 0. };
	} else {
		out = d65_rel_intensity(wavelength) * xyz_from_wavelength(wavelength);
	}
}

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

void print_color(const char * text, const RGB& color) {
	int r = color.r * 255;
	int g = color.g * 255;
	int b = color.b * 255;
	std::cout << "\x1b[38;2;" << r << ";" << g << ";" << b << "m";
	std::cout << text;
	std::cout << "\x1b[0m";
}

int main() {
	// XYZ output = emitted_color(d65_spectrum, rust_func);
	// std::cout << "(" << output.x << ", " << output.y << ", " << output.z << ")\n";
	// RGB s_out = sRGB_from_XYZ(output);
	// std::cout << "(" << s_out.r * 255 << ", " << s_out.g * 255 << ", " << s_out.b * 255 << ")\n";
	std::string out_text;
	out_text.reserve(128);
	for (double i = 1.60; i < 3.5; i += 0.006125) {
		out_text.clear();
		Simple_Gap_Function sg{i};
		XYZ output = emitted_color(d65_spectrum, sg);
		RGB s_out = sRGB_from_XYZ(output);
		out_text += std::to_string(i);
		out_text += " eV";
		print_color(out_text.c_str(), s_out);
		std::cout << "\n";
		// std::cout
		// 	<< i
		// 	<< " eV:\t\t("
		// 	<< s_out.r * 255
		// 	<< ", "
		// 	<< s_out.g * 255
		// 	<< ", "
		// 	<< s_out.b * 255
		// 	<< ")\n\t\t\t\t("
		// 	<< output.x
		// 	<< ", "
		// 	<< output.y
		// 	<< ", "
		// 	<< output.z
		// 	<< ")\n";
	}
	return 0;
}
