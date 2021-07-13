#include "xorshift.h"

uint32_t xorshift32(struct xorshift32_state *state)
{
	/* Algorithm "xor" from p. 4 of Marsaglia, "Xorshift RNGs" */
	uint32_t x = state->a;
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 5;
	return state->a = x;
}

// Unbiased random number range
uint32_t bounded_rand(struct xorshift32_state* state, uint32_t range) {
    uint32_t x = xorshift32(state);
    uint64_t m = (uint64_t)(x) * (uint64_t)(range);
    uint32_t l = (uint32_t)(m);
    if (l < range) {
        uint32_t t = -range;
        if (t >= range) {
            t -= range;
            if (t >= range)
                t %= range;
        }
        while (l < t) {
            x = xorshift32(state);
            m = (uint64_t)(x) * (uint64_t)(range);
            l = (uint32_t)(m);
        }
    }
    return m >> 32;
}


float rand_float(struct xorshift32_state * state)
{
    uint32_t cur_value = xorshift32(state);
    uint32_t max_int = (uint32_t)-1;
    return (float)cur_value / (float)max_int;
};