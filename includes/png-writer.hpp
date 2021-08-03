#pragma once
#include "png-writer.h"
#include <memory>

namespace png {

class PNG {
public:
    PNG(size_t width, size_t height);
    PNG(const PNG& other) = delete;
    PNG& operator=(const PNG& other) = delete;
    PNG(PNG&& other) noexcept;
    PNG& operator=(PNG&& other) noexcept;
    pixel_t& at(size_t row, size_t col);
    pixel_t at(size_t row, size_t col) const;
    int write_to_file(const char * path);
    ~PNG();
private:
    bitmap_t bits;
};

}