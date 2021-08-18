#include "xorshift.hpp"

Xorshift32::Xorshift32(uint32_t seed) : cur_state(seed) {}

uint32_t Xorshift32::get_rand() {
    uint32_t& x = cur_state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    return x;
}

float Xorshift32::rand_float() {
    uint32_t cur_value = get_rand();
    uint32_t max_int = (uint32_t)-1;
    return (float)cur_value / (float)max_int;
}
