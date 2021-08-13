#pragma once
#include "LCh.hpp"
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
    static constexpr double k_h = 1.0;
    static constexpr double k_C = 1.0;
    static constexpr double k_L = 1.0;
};