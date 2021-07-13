#pragma once
#include "XYZ.h"

class Black_Body {
public:
    Black_Body(double T);
    double operator()(double wavelength);
    void operator()(double& out, double wavelength);
    void operator()(XYZ& out, double wavelength);
private:
    double temperature;
};