#include "shortest-path-color.hpp"
#include "color-difference-2000.hpp"
#include "print-color.hpp"
#include <unordered_map>
#include <cmath>
#include <queue>
#include <limits>
#include <string>
#include <algorithm>
#include <iostream>

namespace std {

  template <>
  struct hash<LCh>
  {
    std::size_t operator()(const LCh& lch) const
    {
        int L = round(lch.L * 100.0);
        int C = round(lch.C * 100.0);
        int h = round((180.0 / M_PI) * lch.h * 100.0);
        std::string out;
        out.reserve(64);
        out += std::to_string(L);
        out += " ";
        out += std::to_string(C);
        out += " ";
        out += std::to_string(h);
        std::hash<std::string> hasher;
        return hasher(out);
    }
  };

}

class LCh_Hash {
public:
    std::string operator()(const LCh& lch) {
        int L = round(lch.L * 100.0);
        int C = round(lch.C * 100.0);
        int h = round((180.0 / M_PI) * lch.h * 100.0);
        std::string out;
        out.reserve(64);
        out += std::to_string(L);
        out += " ";
        out += std::to_string(C);
        out += " ";
        out += std::to_string(h);
        return out;
    }
};

struct Color_Node {
    LCh prev_color;
    double dist_estimate;
    double dist_traveled;
	double dist_from_prev;
    bool operator<(const Color_Node& other) const {
        return dist_estimate > other.dist_estimate;
    }
};

static bool operator==(const LCh& a, const LCh& b) {
    static LCh_Hash checker;
    return checker(a) == checker(b);
}

std::vector<LCh> shortest_path_approx(const LCh& start, const LCh& stop) {
    int L1 = round(start.L * 10);
    int L2 = round(stop.L * 10);
    int C1 = round(start.C * 10);
    int C2 = round(stop.C * 10);
    int h1 = round(((180.0 / M_PI) * start.h) * 10);
    int h2 = round(((180.0 / M_PI) * stop.h) * 10);
    std::vector<LCh> out;
    out.reserve(abs(L1 - L2) + abs(C1 - C2) + abs(h1 - h2) + 2);
    LCh_Hash checker;
    std::unordered_map<LCh, Color_Node> visited;
    visited.reserve(20000000);
    auto priority = [&visited](const LCh& a, const LCh& b){ return visited[a].dist_estimate > visited[b].dist_estimate; };
    std::priority_queue<LCh, std::vector<LCh>, decltype(priority)> nodes(priority);
    nodes.emplace(start);
    std::array<std::array<double, 3>, 6> cube{
        std::array<double, 3>{0.1, 0.0, 0.0},
        std::array<double, 3>{0.0, 0.1, 0.0},
        std::array<double, 3>{0.0, 0.0, 0.1 * M_PI / 180.0},
        std::array<double, 3>{-0.1, 0.0, 0.0},
        std::array<double, 3>{0.0, -0.1, 0.0},
        std::array<double, 3>{0.0, 0.0, -0.1 * M_PI / 180.0}
    };
    visited[start] = {
        {0, 0, 0},
        color_distance(start, stop),
        0.0,
		0.0
    };
    std::string stop_hash = checker(stop);
    double min_dist = std::numeric_limits<double>::infinity();
    double scale = 5.0;
    std::cout << "Cutoff Dist: " << 3.0005 * scale * 0.1 << "\n";
	LCh cur_color = nodes.top();
    while (!nodes.empty() && (min_dist > 3.0005 * scale * 0.1) && (checker(nodes.top()) != stop_hash)) {
		cur_color = nodes.top();
        nodes.pop();
        for (auto a : cube) {
            LCh intermediate_color = cur_color;
            for (size_t i = 0; i < 3; i++) {
                intermediate_color[i] += scale * a[i];
            }
            if (
                intermediate_color.L <= 0.0
                || intermediate_color.L >= 200.0
                || intermediate_color.C <= 0.0
                || intermediate_color.C >= 200.0
                || intermediate_color.h <= -4.0 * M_PI
                || intermediate_color.h >= 4.0 * M_PI
                || fabs(intermediate_color.L - start.L) - scale * 0.1 > fabs(stop.L - start.L)
                || fabs(intermediate_color.C - start.C) - scale * 0.1 > fabs(stop.C - start.C)
                || fabs(intermediate_color.h - start.h) - scale * 0.1 > fabs(stop.h - start.h)
            ) {
                continue;
            }
            double dist_step = color_distance(intermediate_color, start);
            double dist_traveled = visited[cur_color].dist_traveled + dist_step;
            double dist_estimate = visited[cur_color].dist_traveled
                + color_distance(intermediate_color, stop);
            if (!visited.count(intermediate_color)) {
                visited[intermediate_color] = {
                    cur_color,
                    dist_estimate,
                    dist_traveled,
					dist_step
                };
                nodes.push(intermediate_color);
            } else if (visited[intermediate_color].dist_estimate > dist_estimate) {
                visited[intermediate_color] = {
                    cur_color,
                    dist_estimate,
                    dist_traveled,
					dist_step
                };
            }
            double possible_min_dist = 0.0;
            std::array<double, 3> scale_factors = {
                1.0,
                1.0,
                180.0 / M_PI
            };
            for (size_t i = 0; i < 3; i++) {
                double dist = intermediate_color[i] - stop[i];
                dist *= scale_factors[i];
                possible_min_dist += dist * dist;
            }
            if (possible_min_dist < min_dist) {
                min_dist = possible_min_dist;
                std::cout << "\r\x1b[KMin Dist: " << min_dist;
                full_print_color(intermediate_color);
				if (min_dist < 3.0005 * scale * 0.1) {
					cur_color = intermediate_color;
					break;
				}
            }
        }
    }
    std::cout << "Reached final point!\n";
    std::cout << "Nodes Visited: " << visited.size() << "\n";
	std::cout << "Distance Traveled: " << visited[cur_color].dist_traveled << "\n";
	std::cout << "Min Possible Dist: " << visited[start].dist_estimate << "\n";
	double min_dist_traveled = std::numeric_limits<double>::infinity();
    while (visited[cur_color].dist_traveled > 0.0001) {
		if (visited[cur_color].dist_from_prev < min_dist_traveled) {
			min_dist_traveled = visited[cur_color].dist_from_prev;
		}
        out.push_back(cur_color);
        cur_color = visited[cur_color].prev_color;
    }
	std::vector<LCh> adj_out;
	adj_out.reserve(out.size() * 2);
	adj_out.push_back(stop);
	for (auto & c : out) {
		size_t ratio = visited[c].dist_from_prev / min_dist_traveled;
		for (size_t i = 0; i <= ratio; i++) {
			adj_out.push_back(close_interpolate(visited[c].prev_color, c, (double)i / (double)ratio));
		}
	}
	adj_out.push_back(start);
    std::reverse(adj_out.begin(), adj_out.end());
	for (size_t i = 0; i < adj_out.size(); i++) {
		const auto& a = adj_out[i];
		if (
			a.C <= 0.0
			|| a.C >= 200.0
			|| a.L <= 0.0
			|| a.L >= 200.0
			|| a.h <= - 6.0 * M_PI
			|| a.h >= 6.0 * M_PI
			|| std::isnan(a.L)
			|| std::isnan(a.C)
			|| std::isnan(a.h)
		) {
			std::cout << "ADJ " << i << ": (" << a.L << ", " << a.C << ", " << a.h << ")\n";
		}
	}
    return adj_out;
}
