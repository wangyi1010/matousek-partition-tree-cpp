#include "matousek_partition_tree/core.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <numeric>
#include <sstream>
#include <utility>

namespace mpt {
namespace {

Rational approximate(double value) {
    constexpr std::int64_t scale = 1'000'000'000;
    return rational(static_cast<std::int64_t>(std::llround(value * scale)), scale);
}

Point segment_line_point(const Point &a, const Point &b, const Line &line) {
    const Rational sa = side(line, a);
    const Rational sb = side(line, b);
    const Rational t = sa / (sa - sb);
    return {a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)};
}

template <class T> bool contains(const std::vector<T> &values, const T &value) {
    return std::find(values.begin(), values.end(), value) != values.end();
}

std::size_t weighted_index(const std::vector<std::size_t> &indices,
                           const std::vector<double> &weights, std::mt19937_64 &rng) {
    std::vector<double> selected;
    selected.reserve(indices.size());
    for (const auto index : indices)
        selected.push_back(weights[index]);
    std::discrete_distribution<std::size_t> distribution(selected.begin(), selected.end());
    return indices[distribution(rng)];
}

} // namespace

Rational rational(std::int64_t numerator, std::int64_t denominator) {
    return Rational(Integer(numerator), Integer(denominator));
}

double to_double(const Rational &value) {
    return value.numerator().convert_to<double>() / value.denominator().convert_to<double>();
}

Line dual_of_point(const Point &point) { return {point.x, -point.y}; }
Point dual_of_line(const Line &line) { return {line.m, -line.c}; }

Rational side(const Line &line, const Point &point) {
    return point.y - (line.m * point.x + line.c);
}

Polygon clip(const Polygon &polygon, const Line &line, bool keep_above) {
    if (polygon.empty())
        return {};
    Polygon out;
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        const Point &a = polygon[i];
        const Point &b = polygon[(i + 1) % polygon.size()];
        const Rational sa = side(line, a);
        const Rational sb = side(line, b);
        const bool in_a = keep_above ? sa >= 0 : sa <= 0;
        const bool in_b = keep_above ? sb >= 0 : sb <= 0;
        if (in_a)
            out.push_back(a);
        if (in_a != in_b && sa != sb)
            out.push_back(segment_line_point(a, b, line));
    }
    Polygon deduplicated;
    for (std::size_t i = 0; i < out.size(); ++i) {
        if (out[i] != out[(i + out.size() - 1) % out.size()])
            deduplicated.push_back(out[i]);
    }
    return deduplicated.size() >= 3 ? deduplicated : Polygon{};
}

Rational poly_area2(const Polygon &polygon) {
    Rational sum = rational(0);
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        const Point &a = polygon[i];
        const Point &b = polygon[(i + 1) % polygon.size()];
        sum += a.x * b.y - b.x * a.y;
    }
    return sum < 0 ? -sum : sum;
}

std::vector<Polygon> arrangement_cells(const std::vector<Line> &lines, const Polygon &box) {
    std::vector<Polygon> cells{box};
    for (const auto &line : lines) {
        std::vector<Polygon> next;
        for (const auto &cell : cells) {
            for (bool above : {true, false}) {
                auto part = clip(cell, line, above);
                if (!part.empty() && poly_area2(part) > 0)
                    next.push_back(std::move(part));
            }
        }
        cells = std::move(next);
    }
    return cells;
}

std::vector<Triangle> fan_triangles(const Polygon &polygon) {
    std::vector<Triangle> triangles;
    for (std::size_t i = 1; i + 1 < polygon.size(); ++i) {
        Triangle triangle{polygon[0], polygon[i], polygon[i + 1]};
        if (poly_area2(triangle) > 0)
            triangles.push_back(std::move(triangle));
    }
    return triangles;
}

bool line_crosses_tri(const Line &line, const Triangle &triangle) {
    bool negative = false;
    bool positive = false;
    for (const auto &vertex : triangle) {
        const auto value = side(line, vertex);
        negative = negative || value < 0;
        positive = positive || value > 0;
    }
    return negative && positive;
}

bool point_in_tri(const Point &point, const Triangle &triangle) {
    assert(triangle.size() == 3);
    const auto cross = [&](const Point &a, const Point &b) {
        return (b.x - a.x) * (point.y - a.y) - (b.y - a.y) * (point.x - a.x);
    };
    const Rational d1 = cross(triangle[0], triangle[1]);
    const Rational d2 = cross(triangle[1], triangle[2]);
    const Rational d3 = cross(triangle[2], triangle[0]);
    const bool negative = d1 < 0 || d2 < 0 || d3 < 0;
    const bool positive = d1 > 0 || d2 > 0 || d3 > 0;
    return !(negative && positive);
}

Polygon combinatorial_box(const std::vector<Line> &lines, const std::vector<Point> &points) {
    bool initialized = false;
    Rational min_x = rational(0);
    Rational max_x = rational(0);
    Rational min_y = rational(0);
    Rational max_y = rational(0);
    const auto include = [&](const Rational &x, const Rational &y) {
        if (!initialized) {
            min_x = max_x = x;
            min_y = max_y = y;
            initialized = true;
            return;
        }
        min_x = std::min(min_x, x);
        max_x = std::max(max_x, x);
        min_y = std::min(min_y, y);
        max_y = std::max(max_y, y);
    };
    for (const auto &point : points) {
        include(point.x, point.y);
    }
    if (!lines.empty()) {
        // For every nonparallel pair, the numerator of the intersection's
        // x-coordinate is at most the intercept span and its denominator is at
        // least the smallest nonzero slope gap. This exact conservative bound
        // contains every arrangement vertex without enumerating all pairs.
        const auto absolute = [](const Rational &value) { return value < 0 ? -value : value; };
        Rational min_c = lines.front().c;
        Rational max_c = lines.front().c;
        Rational max_abs_m = absolute(lines.front().m);
        Rational max_abs_c = absolute(lines.front().c);
        std::vector<Rational> slopes;
        slopes.reserve(lines.size());
        for (const auto &line : lines) {
            slopes.push_back(line.m);
            min_c = std::min(min_c, line.c);
            max_c = std::max(max_c, line.c);
            max_abs_m = std::max(max_abs_m, absolute(line.m));
            max_abs_c = std::max(max_abs_c, absolute(line.c));
        }
        std::sort(slopes.begin(), slopes.end());
        std::optional<Rational> minimum_slope_gap;
        for (std::size_t i = 1; i < slopes.size(); ++i) {
            const Rational gap = slopes[i] - slopes[i - 1];
            if (gap > 0 && (!minimum_slope_gap || gap < *minimum_slope_gap))
                minimum_slope_gap = gap;
        }
        const Rational x_bound =
            minimum_slope_gap ? (max_c - min_c) / *minimum_slope_gap : rational(0);
        const Rational y_bound = max_abs_m * x_bound + max_abs_c;
        include(-x_bound, -y_bound);
        include(x_bound, y_bound);
    }
    if (!initialized)
        include(rational(0), rational(0));
    const Rational cx = (min_x + max_x) / rational(2);
    const Rational cy = (min_y + max_y) / rational(2);
    const Rational half = std::max({max_x - min_x, max_y - min_y, rational(1)}) * rational(2);
    const Rational lo_x = cx - half;
    const Rational hi_x = cx + half;
    const Rational lo_y = cy - half;
    const Rational hi_y = cy + half;
    return {{lo_x, lo_y}, {hi_x, lo_y}, {hi_x, hi_y}, {lo_x, hi_y}};
}

bool poly_crossed(const Line &line, const Polygon &polygon) {
    bool negative = false;
    bool positive = false;
    for (const auto &vertex : polygon) {
        const auto value = side(line, vertex);
        negative = negative || value < 0;
        positive = positive || value > 0;
    }
    return negative && positive;
}

std::vector<Triangle> weighted_cutting(const std::vector<Line> &lines,
                                       const std::vector<Rational> &weights, double t,
                                       const Polygon &box, std::mt19937_64 &rng,
                                       std::optional<std::size_t> max_faces, int max_tries) {
    if (lines.size() != weights.size())
        throw std::invalid_argument("lines/weights size mismatch");
    const Rational total = std::accumulate(weights.begin(), weights.end(), rational(0));
    const Rational t_fraction = approximate(t);
    if (lines.empty() || t <= 1.0) {
        auto triangles = fan_triangles(box);
        if (max_faces && triangles.size() > *max_faces) {
            throw CuttingError("trivial cutting exceeds face budget");
        }
        return triangles;
    }
    std::vector<double> floating_weights;
    for (const auto &weight : weights)
        floating_weights.push_back(to_double(weight));

    auto heavy = [&](const Polygon &cell) {
        std::vector<std::size_t> indices;
        Rational crossing_weight = rational(0);
        for (std::size_t i = 0; i < lines.size(); ++i) {
            if (poly_crossed(lines[i], cell)) {
                indices.push_back(i);
                crossing_weight += weights[i];
            }
        }
        return std::pair{crossing_weight * t_fraction > total, indices};
    };

    std::size_t k0 = std::max<std::size_t>(1, std::min<std::size_t>(lines.size(), std::ceil(t)));
    for (int attempt = 0; attempt < max_tries; ++attempt) {
        std::vector<Line> sample;
        if (k0 > 0) {
            std::vector<std::size_t> all(lines.size());
            std::iota(all.begin(), all.end(), 0);
            for (std::size_t draw = 0; draw < 4 * k0 && sample.size() < k0; ++draw) {
                const auto index = weighted_index(all, floating_weights, rng);
                if (!contains(sample, lines[index]))
                    sample.push_back(lines[index]);
            }
        }
        std::vector<Polygon> stack = arrangement_cells(sample, box);
        std::vector<Polygon> final_cells;
        bool failed = false;
        int guard = 0;
        while (!stack.empty()) {
            if (++guard > 4000) {
                failed = true;
                break;
            }
            Polygon cell = std::move(stack.back());
            stack.pop_back();
            auto [is_heavy, indices] = heavy(cell);
            if (!is_heavy) {
                final_cells.push_back(std::move(cell));
                continue;
            }
            const std::size_t splitter = weighted_index(indices, floating_weights, rng);
            auto above = clip(cell, lines[splitter], true);
            auto below = clip(cell, lines[splitter], false);
            int pieces = 0;
            for (auto *part : {&above, &below}) {
                if (!part->empty() && poly_area2(*part) > 0) {
                    stack.push_back(std::move(*part));
                    ++pieces;
                }
            }
            if (pieces < 2) {
                failed = true;
                break;
            }
        }
        if (failed)
            continue;
        std::vector<Triangle> triangles;
        for (const auto &cell : final_cells) {
            auto pieces = fan_triangles(cell);
            triangles.insert(triangles.end(), pieces.begin(), pieces.end());
        }
        if (max_faces && triangles.size() > *max_faces) {
            k0 = k0 > 0 ? k0 - std::max<std::size_t>(1, k0 / 3) : 0;
            continue;
        }
        return triangles;
    }
    std::ostringstream message;
    message << "no verified (1/" << t << ")-cutting in " << max_tries << " tries";
    throw CuttingError(message.str());
}

std::vector<Line> build_test_set(const std::vector<Point> &points, double r, std::mt19937_64 &rng) {
    std::vector<Line> duals;
    for (const auto &point : points)
        duals.push_back(dual_of_point(point));
    const auto triangles = weighted_cutting(duals, std::vector<Rational>(duals.size(), rational(1)),
                                            0.25 * std::sqrt(r), combinatorial_box(duals), rng);
    std::vector<Point> vertices;
    for (const auto &triangle : triangles) {
        for (const auto &vertex : triangle) {
            if (!contains(vertices, vertex))
                vertices.push_back(vertex);
        }
    }
    const auto budget = std::max<std::size_t>(4, std::ceil(r));
    if (vertices.size() > budget) {
        std::ostringstream message;
        message << "cutting has " << vertices.size() << " vertices, exceeding r=" << budget;
        throw TestSetError(message.str());
    }
    std::vector<Line> result;
    for (const auto &vertex : vertices)
        result.push_back(dual_of_point(vertex));
    return result;
}

PartitionResult simplicial_partition(const std::vector<Point> &points, std::size_t s,
                                     std::mt19937_64 &rng, double cutting_scale) {
    if (s == 0)
        throw std::invalid_argument("s must be positive");
    if (!(cutting_scale > 0.0))
        throw std::invalid_argument("cutting_scale must be positive");
    std::vector<Line> Q = build_test_set(points, static_cast<double>(points.size()) / s, rng);
    std::vector<int> kappa(Q.size(), 0);
    std::vector<Point> remaining = points;
    const Polygon box = combinatorial_box(Q, points);
    std::vector<Group> groups;
    while (remaining.size() >= 2 * s) {
        const std::size_t face_budget = remaining.size() / s;
        const double t = cutting_scale * std::sqrt(static_cast<double>(remaining.size()) / s);
        std::vector<Rational> weights;
        for (const int exponent : kappa) {
            Integer value = Integer(1) << exponent;
            weights.emplace_back(value);
        }
        const auto triangles = weighted_cutting(Q, weights, t, box, rng, face_budget);
        std::vector<std::size_t> counts;
        for (const auto &triangle : triangles) {
            counts.push_back(std::count_if(remaining.begin(), remaining.end(), [&](const Point &p) {
                return point_in_tri(p, triangle);
            }));
        }
        const auto best = static_cast<std::size_t>(
            std::distance(counts.begin(), std::max_element(counts.begin(), counts.end())));
        if (counts[best] < s)
            throw std::logic_error("pigeonhole invariant violated");
        const Triangle triangle = triangles[best];
        std::vector<Point> inside;
        std::vector<Point> outside;
        for (const auto &point : remaining) {
            (point_in_tri(point, triangle) ? inside : outside).push_back(point);
        }
        std::vector<Point> group(inside.begin(), inside.begin() + static_cast<std::ptrdiff_t>(s));
        outside.insert(outside.end(), inside.begin() + static_cast<std::ptrdiff_t>(s),
                       inside.end());
        remaining = std::move(outside);
        groups.push_back({std::move(group), triangle});
        for (std::size_t i = 0; i < Q.size(); ++i) {
            if (line_crosses_tri(Q[i], triangle))
                ++kappa[i];
        }
    }
    if (!remaining.empty())
        groups.push_back({remaining, bounding_triangle(remaining)});
    const int maximum = kappa.empty() ? 0 : *std::max_element(kappa.begin(), kappa.end());
    return {std::move(groups), {std::move(Q), std::move(kappa), maximum}};
}

Triangle bounding_triangle(const std::vector<Point> &points) {
    if (points.empty())
        throw std::invalid_argument("cannot bound an empty point set");
    Rational lo_x = points.front().x, hi_x = points.front().x;
    Rational lo_y = points.front().y, hi_y = points.front().y;
    for (const auto &point : points) {
        lo_x = std::min(lo_x, point.x);
        hi_x = std::max(hi_x, point.x);
        lo_y = std::min(lo_y, point.y);
        hi_y = std::max(hi_y, point.y);
    }
    const Rational width = hi_x - lo_x + rational(1);
    const Rational height = hi_y - lo_y + rational(1);
    return {{lo_x - width, lo_y - rational(1)},
            {hi_x + width, lo_y - rational(1)},
            {(lo_x + hi_x) / rational(2), hi_y + rational(2) * height}};
}

std::unique_ptr<PNode> build_tree(const std::vector<Point> &points, std::size_t r,
                                  std::size_t leaf_size, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    return build_tree(points, r, leaf_size, rng);
}

std::unique_ptr<PNode> build_tree(const std::vector<Point> &points, std::size_t r,
                                  std::size_t leaf_size, std::mt19937_64 &rng,
                                  std::optional<Triangle> simplex) {
    if (r < 2 || leaf_size == 0)
        throw std::invalid_argument("r >= 2 and leaf_size > 0 required");
    auto node = std::make_unique<PNode>();
    node->count = points.size();
    node->simplex = std::move(simplex);
    if (points.size() <= leaf_size) {
        node->points = points;
        return node;
    }
    const std::size_t s = std::max<std::size_t>(1, points.size() / r);
    auto partition = simplicial_partition(points, s, rng);
    for (const auto &group : partition.groups) {
        node->children.push_back(build_tree(group.points, r, leaf_size, rng, group.simplex));
    }
    return node;
}

Rational halfplane_side(const Halfplane &halfplane, const Point &point) {
    return halfplane.a * point.x + halfplane.b * point.y + halfplane.c;
}

std::size_t query_count(const PNode &node, const Halfplane &halfplane) {
    if (node.simplex) {
        bool all_inside = true;
        bool all_outside = true;
        for (const auto &vertex : *node.simplex) {
            const auto value = halfplane_side(halfplane, vertex);
            all_inside = all_inside && value >= 0;
            all_outside = all_outside && value < 0;
        }
        if (all_inside)
            return node.count;
        if (all_outside)
            return 0;
    }
    if (node.points) {
        return std::count_if(node.points->begin(), node.points->end(), [&](const Point &point) {
            return halfplane_side(halfplane, point) >= 0;
        });
    }
    std::size_t result = 0;
    for (const auto &child : node.children)
        result += query_count(*child, halfplane);
    return result;
}

std::size_t crossing_number(const PNode &node, const Line &line) {
    return std::count_if(node.children.begin(), node.children.end(), [&](const auto &child) {
        return child->simplex && line_crosses_tri(line, *child->simplex);
    });
}

} // namespace mpt
