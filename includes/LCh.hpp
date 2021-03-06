#pragma once
#include "XYZ.hpp"

struct LCh {
public:
    LCh();
    LCh(double luma, double chroma, double hue);
    LCh(XYZ in);
    XYZ to_XYZ() const;
    double& operator[](size_t i);
    double operator[](size_t i) const;
public:
    double L, C, h;
};

LCh close_interpolate(const LCh& a, const LCh& b, double t);
LCh interpolate(const LCh& a, const LCh& b, double t);
