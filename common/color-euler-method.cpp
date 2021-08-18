#include "color-euler-method.hpp"
#include "color-christoffel-symbols.hpp"
#include <iostream>
#include <vector>
#include <cstddef>
#include <cmath>
#include <array>

class RunningStat {
public:
    RunningStat() : m_n(0) {}
    void Clear() {
        m_n = 0;
    }
    void Push(double x) {
        m_n++;
        // See Knuth TAOCP vol 2, 3rd edition, page 232
        if (m_n == 1) {
            m_oldM = m_newM = x;
            m_oldS = 0.0;
        } else {
            m_newM = m_oldM + (x - m_oldM)/m_n;
            m_newS = m_oldS + (x - m_oldM)*(x - m_newM);

            // set up for next iteration
            m_oldM = m_newM; 
            m_oldS = m_newS;
        }
    }
    int NumDataValues() const {
        return m_n;
    }
    double Mean() const {
        return (m_n > 0) ? m_newM : 0.0;
    }
    double Variance() const {
        return ( (m_n > 1) ? m_newS/(m_n - 1) : 0.0 );
    }
    double StandardDeviation() const {
        return sqrt( Variance() );
    }
private:
    int m_n;
    double m_oldM, m_newM, m_oldS, m_newS;
};

std::vector<LCh> out_color(LCh start, double dL_0, double dC_0, double dh_0, double step, size_t num_elements) {
    RunningStat tracker_L, tracker_C, tracker_h;
    Color_Christoffel_Symbols gamma;
    std::vector<LCh> out;
    out.reserve(num_elements);
    out.push_back(start);
    LCh pos = start;
    std::array<double, 3> acc{0, 0, 0};
    double ddL, ddC, ddh;
    double norm = dL_0 * dL_0
        + dC_0 * dC_0
        + dh_0 * dh_0;
    norm = sqrt(norm);
    std::array<double, 3> vel{
        dL_0 / norm,
        dC_0 / norm,
        dh_0 / norm
    };
    for (size_t n = 1; n < num_elements; n++) {
        gamma.calculate(pos);
        acc[0] = acc[1] = acc[2] = 0.0;
        for (size_t i = 0; i < 3; i++) {
            for (size_t j = 0; j < 3; j++) {
                for (size_t k = 0; k < 3; k++) {
                    acc[i] -= gamma.gamma[i][j][k] * vel[j] * vel[k];
                }
            }
        }
        for (size_t i = 0; i < 3; i++) {
            vel[i] += acc[i] * step;
        }
        norm = 0.0;
        for (const auto& v : vel) {
            norm += v * v;
        }
        norm = sqrt(norm);
        for (auto& v : vel) {
            v /= norm;
        }
        tracker_L.Push(vel[0]);
        tracker_C.Push(vel[1]);
        tracker_h.Push(vel[2]);
        for (size_t i = 0; i < 3; i++) {
            pos[i] += vel[i] * step;
        }
        out.push_back(pos);
    }
    return out;
}