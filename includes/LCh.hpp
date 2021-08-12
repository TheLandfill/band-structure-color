#pragma once
#include "XYZ.hpp"

struct LCh {
public:
    LCh(double luma, double chroma, double hue);
    XYZ to_XYZ();
public:
    double L, C, h;
};