#include "matousek_partition_tree/core.hpp"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <vector>

namespace {
int failures = 0;

void check(bool condition, const std::string &message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

std::vector<mpt::Point> make_points(std::size_t n = 400, std::uint64_t seed = 42) {
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<std::int64_t> distribution(0, 10'000);
    std::vector<mpt::Point> points;
    for (std::size_t i = 0; i < n; ++i) {
        points.push_back(
            {mpt::rational(distribution(rng), 10'000), mpt::rational(distribution(rng), 10'000)});
    }
    return points;
}

void test_partition() {
    auto points = make_points();
    std::mt19937_64 rng(42);
    auto result = mpt::simplicial_partition(points, 16, rng);
    check(!result.stats.Q.empty(), "nontrivial partition builds a test set");
    std::vector<mpt::Point> seen;
    for (const auto &group : result.groups) {
        check(group.points.size() >= 16 && group.points.size() < 32, "group size is in [s, 2s)");
        for (const auto &point : group.points) {
            check(mpt::point_in_tri(point, group.simplex), "group point lies in its simplex");
            seen.push_back(point);
        }
    }
    check(seen.size() == points.size(), "partition covers every point");
    for (const auto &point : points) {
        check(std::count(seen.begin(), seen.end(), point) ==
                  std::count(points.begin(), points.end(), point),
              "partition preserves point multiplicity");
    }
    for (std::size_t i = 0; i < result.stats.Q.size(); ++i) {
        int recount = 0;
        for (std::size_t j = 0; j + 1 < result.groups.size(); ++j) {
            recount += mpt::line_crosses_tri(result.stats.Q[i], result.groups[j].simplex);
        }
        check(recount == result.stats.kappa[i], "test-set crossing bookkeeping is exact");
    }
}

void test_cutting() {
    std::mt19937_64 rng(3);
    std::uniform_int_distribution<std::int64_t> distribution(-100, 100);
    std::vector<mpt::Line> lines;
    std::vector<mpt::Rational> weights;
    for (int i = 0; i < 30; ++i) {
        lines.push_back(
            {mpt::rational(distribution(rng), 100), mpt::rational(distribution(rng), 100)});
        weights.push_back(mpt::rational(std::int64_t{1} << (i % 6)));
    }
    mpt::Polygon box{{mpt::rational(-4), mpt::rational(-4)},
                     {mpt::rational(4), mpt::rational(-4)},
                     {mpt::rational(4), mpt::rational(4)},
                     {mpt::rational(-4), mpt::rational(4)}};
    auto triangles = mpt::weighted_cutting(lines, weights, 3.0, box, rng, 40);
    check(triangles.size() <= 40, "cutting respects face budget");
    const auto total = std::accumulate(weights.begin(), weights.end(), mpt::rational(0));
    for (const auto &triangle : triangles) {
        auto crossing = mpt::rational(0);
        for (std::size_t i = 0; i < lines.size(); ++i) {
            if (mpt::line_crosses_tri(lines[i], triangle))
                crossing += weights[i];
        }
        check(crossing * mpt::rational(3) <= total, "cutting cell is light");
    }
    bool threw = false;
    try {
        (void)mpt::weighted_cutting({}, {}, 1.0, box, rng, 1);
    } catch (const mpt::CuttingError &) {
        threw = true;
    }
    check(threw, "impossible face budget throws CuttingError");
}

void test_combinatorial_box() {
    const std::vector<mpt::Line> lines{
        {mpt::rational(-3, 2), mpt::rational(7, 3)}, {mpt::rational(-1, 4), mpt::rational(-5, 2)},
        {mpt::rational(2, 3), mpt::rational(11, 5)}, {mpt::rational(2, 3), mpt::rational(-9, 7)},
        {mpt::rational(5, 2), mpt::rational(1, 9)},
    };
    const std::vector<mpt::Point> points{{mpt::rational(-20), mpt::rational(4)},
                                         {mpt::rational(3), mpt::rational(30)}};
    const auto box = mpt::combinatorial_box(lines, points);
    const auto lo_x = box[0].x;
    const auto lo_y = box[0].y;
    const auto hi_x = box[2].x;
    const auto hi_y = box[2].y;
    for (const auto &point : points) {
        check(lo_x <= point.x && point.x <= hi_x && lo_y <= point.y && point.y <= hi_y,
              "combinatorial box contains input points");
    }
    for (std::size_t i = 0; i < lines.size(); ++i) {
        for (std::size_t j = i + 1; j < lines.size(); ++j) {
            if (lines[i].m == lines[j].m)
                continue;
            const auto x = (lines[j].c - lines[i].c) / (lines[i].m - lines[j].m);
            const auto y = lines[i].m * x + lines[i].c;
            check(lo_x <= x && x <= hi_x && lo_y <= y && y <= hi_y,
                  "combinatorial box contains every arrangement vertex");
        }
    }
}

void test_invalid_cutting_scale() {
    auto points = make_points(20, 5);
    std::mt19937_64 rng(5);
    bool threw = false;
    try {
        (void)mpt::simplicial_partition(points, 4, rng, 0.0);
    } catch (const std::invalid_argument &) {
        threw = true;
    }
    check(threw, "nonpositive cutting scale is rejected");
}

void test_duality_and_queries() {
    const std::vector<mpt::Point> values{{mpt::rational(2, 3), mpt::rational(-5, 7)},
                                         {mpt::rational(-11, 13), mpt::rational(17, 19)},
                                         {mpt::rational(0), mpt::rational(0)}};
    for (const auto &point : values) {
        check(mpt::dual_of_line(mpt::dual_of_point(point)) == point, "duality is an involution");
    }
    const std::vector<mpt::Line> lines{{mpt::rational(3, 5), mpt::rational(-7, 11)},
                                       {mpt::rational(-2), mpt::rational(9, 4)},
                                       {mpt::rational(0), mpt::rational(0)}};
    for (const auto &line : lines) {
        check(mpt::dual_of_point(mpt::dual_of_line(line)) == line, "line duality is an involution");
    }

    const std::vector<mpt::Point> boundary_points{{mpt::rational(-2), mpt::rational(0)},
                                                  {mpt::rational(-1), mpt::rational(1)},
                                                  {mpt::rational(0), mpt::rational(0)},
                                                  {mpt::rational(1), mpt::rational(-1)},
                                                  {mpt::rational(2), mpt::rational(0)}};
    auto boundary_tree = mpt::build_tree(boundary_points, 4, 2, 11);
    const std::vector<mpt::Halfplane> boundary_queries{
        {mpt::rational(1), mpt::rational(0), mpt::rational(0)},
        {mpt::rational(0), mpt::rational(1), mpt::rational(0)},
        {mpt::rational(-1), mpt::rational(0), mpt::rational(1)},
        {mpt::rational(1), mpt::rational(1), mpt::rational(0)}};
    for (const auto &query : boundary_queries) {
        const auto expected =
            std::count_if(boundary_points.begin(), boundary_points.end(), [&](const auto &point) {
                return mpt::halfplane_side(query, point) >= 0;
            });
        check(mpt::query_count(*boundary_tree, query) == static_cast<std::size_t>(expected),
              "axis-aligned and boundary queries equal brute force");
    }
    for (const auto seed : {0ULL, 7ULL, 19ULL}) {
        const auto points = make_points(120, seed);
        auto tree = mpt::build_tree(points, 25, 16, seed);
        std::mt19937_64 rng(seed + 100);
        std::uniform_int_distribution<std::int64_t> distribution(-20, 20);
        for (int i = 0; i < 8; ++i) {
            mpt::Halfplane query{mpt::rational(distribution(rng)), mpt::rational(distribution(rng)),
                                 mpt::rational(distribution(rng))};
            if (query.a == 0 && query.b == 0)
                query.a = mpt::rational(1);
            const auto expected =
                std::count_if(points.begin(), points.end(), [&](const auto &point) {
                    return mpt::halfplane_side(query, point) >= 0;
                });
            check(mpt::query_count(*tree, query) == static_cast<std::size_t>(expected),
                  "tree query equals brute force");
        }
    }
}
} // namespace

int main() {
    try {
        test_partition();
        test_cutting();
        test_combinatorial_box();
        test_invalid_cutting_scale();
        test_duality_and_queries();
    } catch (const std::exception &error) {
        std::cerr << "UNCAUGHT: " << error.what() << '\n';
        return 2;
    }
    if (failures == 0)
        std::cout << "All property tests passed.\n";
    return failures == 0 ? 0 : 1;
}
