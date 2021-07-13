#ifndef INTEGRATION_H
#define INTEGRATION_H
#include <cstddef>

struct Simpson_Points {
public:
	double start;
	double dx;
	const size_t n;
	Simpson_Points(double a, double b, size_t n);
};

template <typename T>
struct Simpson_Sums {
public:
    Simpson_Sums() : trapezoid_sum(), midpoint_sum(), value_holder() {}
public:
	T trapezoid_sum;
	T midpoint_sum;
	T value_holder;
};

template <typename T, typename FuncType>
T simpson_rule(FuncType f, const Simpson_Points & sp, Simpson_Sums<T>& ss) {
	double offset = 0.0;
	for (size_t i = 0; i < sp.n; i++) {
		f(ss.value_holder, sp.start + offset);
		ss.trapezoid_sum += ss.value_holder;

		f(ss.value_holder, sp.start + offset + sp.dx * 0.5);
		ss.midpoint_sum += ss.value_holder;

		f(ss.value_holder, sp.start + offset + sp.dx);
		ss.trapezoid_sum += ss.value_holder;

		offset += sp.dx;
	}
	ss.trapezoid_sum *= 0.5;
	const double ONE_THIRD = 1.0 / 3.0;
	return (2.0 * ss.midpoint_sum + ss.trapezoid_sum) * ONE_THIRD * sp.dx;
}

#endif