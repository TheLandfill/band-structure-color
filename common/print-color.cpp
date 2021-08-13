#include "print-color.hpp"
#include <iostream>

void print_color(const char * text, const RGB& color) {
	int r = color.r * 255;
	int g = color.g * 255;
	int b = color.b * 255;
	std::cout << "\x1b[38;2;" << r << ";" << g << ";" << b << "m";
	std::cout << text;
	std::cout << "\x1b[0m";
}