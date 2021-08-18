#include "color-fixed-point-finite-method.hpp"
#include "color-christoffel-symbols.hpp"
#include "shortest-path-color.hpp"
#include "xorshift.hpp"
#include <iostream>
#include <cmath>

static void color_iter(std::vector<LCh>& y, std::vector<LCh>& x) {
    Color_Christoffel_Symbols gamma;
    for (size_t iter = 1; iter < x.size() - 1; iter++) {
        gamma.calculate(x[iter]);
        for (size_t i = 0; i < 3; i++) {
            y[iter][i] = (x[iter - 1][i] + x[iter + 1][i]) / 2.0;
            double gamma_terms = 0.0;
            for (size_t j = 0; j < 3; j++) {
                for (size_t k = 0; k < 3; k++) {
                    gamma_terms +=
                        gamma.gamma[i][j][k]
                        * (x[iter + 1][j] - x[iter - 1][j])
                        * (x[iter + 1][k] - x[iter - 1][k]);
                }
            }
            y[iter][i] += gamma_terms / 4.0;
        }
    }
}

static double max_diff(const std::vector<LCh>& x, const std::vector<LCh>& y) {
    double out = 0.0;
    for (size_t i = 0; i < x.size(); i++) {
        double dist = 0.0;
        for (size_t j = 0; j < 3; j++) {
            double diff = (x[i][j] - y[i][j]);
            dist += diff * diff;
        }
        if (dist > out) {
            out = dist;
        }
    }
    return out;
}

static double total_change(const std::vector<LCh>& x, const std::vector<LCh>& y) {
    double out = 0.0;
    for (size_t i = 0; i < x.size(); i++) {
        double dist = 0.0;
        for (size_t j = 0; j < 3; j++) {
            double diff = (x[i][j] - y[i][j]);
            dist += diff * diff;
        }
        out += dist;
    }
    return out;
}

std::vector<LCh> color_finite_element_fixed_point(std::vector<LCh> path_guess) {
    auto out = path_guess;
	//return out;
    std::cout << "Found shortest path!\n";
    std::vector<LCh> intermediate{out};
    size_t num_iterations = 0;
    do {
        color_iter(intermediate, out);
		out = intermediate;
        num_iterations++;
    } while(max_diff(intermediate, out) > 0.0000001 && num_iterations < 1000000);
    if (num_iterations >= 1000000) {
        std::cerr << "WARNING: DID NOT CONVERGE\n";
    }
    std::cout << "Num Iterations: " << num_iterations << "\nMax diff: " << max_diff(intermediate, out) << "\nTotal Change: " << total_change(out, intermediate) << "\n";
    return out;
}
