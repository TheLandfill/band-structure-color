#include "png-writer.hpp"

extern "C" {
    pixel_t * pixel_at(const bitmap_t * bitmap, int x, int y);
    int save_png_to_file(bitmap_t *bitmap, const char *path);
}

namespace png {

PNG::PNG(size_t width, size_t height) {
    bits.width = width;
    bits.height = height;
    bits.pixels = new pixel_t[width * height];
}

PNG::PNG(PNG&& other) noexcept : bits(std::exchange(other.bits, {nullptr, 0, 0})) {}

PNG& PNG::operator=(PNG&& other) noexcept {
    std::swap(bits, other.bits);
    return *this;
}

pixel_t& PNG::at(size_t row, size_t col) {
    return *pixel_at(&bits, col, row);
}

pixel_t PNG::at(size_t row, size_t col) const {
    return *pixel_at(&bits, col, row);
}

int PNG::write_to_file(const char * path) {
    return save_png_to_file(&bits, path);
}

PNG::~PNG() {
    delete[] bits.pixels;
}

}