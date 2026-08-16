#pragma once

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <random>
#include <stdexcept>
#include <vector>

namespace mpt {

using Integer = boost::multiprecision::cpp_int;
using Rational = boost::rational<Integer>;

struct Point {
    Rational x;
    Rational y;
    bool operator==(const Point &) const = default;
};

struct Line {
    Rational m;
    Rational c;
    bool operator==(const Line &) const = default;
};

struct Halfplane {
    Rational a;
    Rational b;
    Rational c;
};

using Polygon = std::vector<Point>;
using Triangle = Polygon;

class CuttingError : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

class TestSetError : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

struct PartitionStats {
    std::vector<Line> Q;
    std::vector<int> kappa;
    int K_Q{};
};

struct Group {
    std::vector<Point> points;
    Triangle simplex;
};

struct PartitionResult {
    std::vector<Group> groups;
    PartitionStats stats;
};

struct PNode {
    std::size_t count{};
    std::optional<Triangle> simplex;
    std::vector<std::unique_ptr<PNode>> children;
    std::optional<std::vector<Point>> points;
};

Rational rational(std::int64_t numerator, std::int64_t denominator = 1);
double to_double(const Rational &value);
Line dual_of_point(const Point &point);
Point dual_of_line(const Line &line);
Rational side(const Line &line, const Point &point);
Polygon clip(const Polygon &polygon, const Line &line, bool keep_above);
Rational poly_area2(const Polygon &polygon);
std::vector<Polygon> arrangement_cells(const std::vector<Line> &lines, const Polygon &box);
std::vector<Triangle> fan_triangles(const Polygon &polygon);
bool line_crosses_tri(const Line &line, const Triangle &triangle);
bool point_in_tri(const Point &point, const Triangle &triangle);
Polygon combinatorial_box(const std::vector<Line> &lines, const std::vector<Point> &points = {});
bool poly_crossed(const Line &line, const Polygon &polygon);

std::vector<Triangle> weighted_cutting(const std::vector<Line> &lines,
                                       const std::vector<Rational> &weights, double t,
                                       const Polygon &box, std::mt19937_64 &rng,
                                       std::optional<std::size_t> max_faces = std::nullopt,
                                       int max_tries = 60);

std::vector<Line> build_test_set(const std::vector<Point> &points, double r, std::mt19937_64 &rng);
// cutting_scale is the fixed positive constant in t_i = cutting_scale * sqrt(n_i / s).
PartitionResult simplicial_partition(const std::vector<Point> &points, std::size_t s,
                                     std::mt19937_64 &rng, double cutting_scale = 0.35);
Triangle bounding_triangle(const std::vector<Point> &points);
std::unique_ptr<PNode> build_tree(const std::vector<Point> &points, std::size_t r = 8,
                                  std::size_t leaf_size = 32, std::uint64_t seed = 0);
std::unique_ptr<PNode> build_tree(const std::vector<Point> &points, std::size_t r,
                                  std::size_t leaf_size, std::mt19937_64 &rng,
                                  std::optional<Triangle> simplex = std::nullopt);
Rational halfplane_side(const Halfplane &halfplane, const Point &point);
std::size_t query_count(const PNode &node, const Halfplane &halfplane);
std::size_t crossing_number(const PNode &node, const Line &line);

} // namespace mpt
