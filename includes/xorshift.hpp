#pragma once
#include <cstdint>

class Xorshift32 {
public:
    Xorshift32(uint32_t seed);
    uint32_t get_rand();
    float rand_float();
private:
    uint32_t cur_state;
};
