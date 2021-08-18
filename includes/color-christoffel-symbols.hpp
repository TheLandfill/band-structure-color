#pragma once
#include "LCh.hpp"
#include "color-difference-2000.hpp"
#include <array>

struct Color_Christoffel_Symbols {
public:
    void calculate(const LCh& lch);
public:
    std::array<
        std::array<
            std::array<
                double,
                3
            >,
            3
        >,
        3
    > gamma;
};