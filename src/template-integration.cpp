#include "template-integration.hpp"

Simpson_Points::Simpson_Points(double a, double b, size_t n) :
start(a),
dx((b - a)/((double)n)),
n(n) {}
