#pragma once
#include <stdint.h>
#include <stddef.h>

/* A coloured pixel. */

typedef struct {
    uint8_t red;
    uint8_t green;
    uint8_t blue;
    uint8_t alpha;
} pixel_t;

/* A picture. */

typedef struct {
    pixel_t *pixels;
    size_t width;
    size_t height;
} bitmap_t;