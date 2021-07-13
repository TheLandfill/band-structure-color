#pragma once
#include <stdint.h>

struct xorshift32_state {
  uint32_t a;
};

uint32_t xorshift32(struct xorshift32_state* state);

uint32_t bounded_rand(struct xorshift32_state* state, uint32_t range);

float rand_float(struct xorshift32_state * state);

#define SHUFFLE_FUNC_DECL(X) void shuffle_array_##X(X * array, int size, struct xorshift32_state* state);
#define SHUFFLE_FUNC_DEFN(X) void shuffle_array_##X(X * array, int size, struct xorshift32_state* state) { \
	X temp; \
	for (int i = size - 1; i > 0; i--) { \
		int32_t swap_index = bounded_rand(state, i + 1); \
		temp = array[i]; \
		array[i] = array[swap_index]; \
		array[swap_index] = temp; \
	} \
} \

