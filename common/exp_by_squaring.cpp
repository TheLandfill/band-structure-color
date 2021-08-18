#include "exp-by-squaring.hpp"

double exp_by_squaring(double x, int n) {
    double x_prod = 1.0;
    while (n != 0) {
        if (n & 1) {
            x_prod *= x;
        }
        n >>= 1;
        x *= x;
    }
    return x_prod;
}