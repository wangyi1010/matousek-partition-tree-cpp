#include "matousek_partition_tree/core.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <string>
#include <vector>

namespace {
std::vector<mpt::Point> make_points(std::size_t count, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<std::int64_t> distribution(0, 10'000);
    std::vector<mpt::Point> points;
    for (std::size_t i = 0; i < count; ++i) {
        points.push_back(
            {mpt::rational(distribution(rng), 10'000), mpt::rational(distribution(rng), 10'000)});
    }
    return points;
}
} // namespace

int main(int argc, char **argv) try {
    const std::size_t n = argc > 1 ? std::stoull(argv[1]) : 1200;
    const std::uint64_t seed = argc > 2 ? std::stoull(argv[2]) : 42;
    auto points = make_points(n, seed);
    std::unique_ptr<mpt::PNode> tree;
    std::uint64_t construction_seed = seed;
    std::string last_error;
    for (std::uint64_t offset = 0; offset < 64 && !tree; ++offset) {
        construction_seed = seed + offset;
        try {
            std::mt19937_64 build_rng(construction_seed);
            tree = mpt::build_tree(points, 64, 32, build_rng);
        } catch (const std::runtime_error &error) {
            last_error = error.what();
        }
    }
    if (!tree)
        throw std::runtime_error("all deterministic construction attempts failed: " + last_error);
    std::mt19937_64 rng(seed + 1);
    std::cout << "built: n=" << n << ", r=64, construction seed=" << construction_seed
              << ", root children=" << tree->children.size() << '\n';

    std::uniform_int_distribution<std::int64_t> query_distribution(-10'000, 10'000);
    std::size_t matches = 0;
    for (int i = 0; i < 60; ++i) {
        mpt::Halfplane query{mpt::rational(query_distribution(rng), 10'000),
                             mpt::rational(query_distribution(rng), 10'000),
                             mpt::rational(query_distribution(rng), 10'000)};
        if (query.a == 0 && query.b == 0)
            query.a = mpt::rational(1);
        const auto exact = std::count_if(points.begin(), points.end(), [&](const auto &point) {
            return mpt::halfplane_side(query, point) >= 0;
        });
        matches += mpt::query_count(*tree, query) == static_cast<std::size_t>(exact);
    }
    std::cout << "query check vs brute force: " << matches << "/60 exact matches\n";

    std::vector<std::size_t> crossings;
    for (int i = 0; i < 300; ++i) {
        mpt::Line line{mpt::rational(query_distribution(rng), 5'000),
                       mpt::rational(query_distribution(rng), 5'000)};
        crossings.push_back(mpt::crossing_number(*tree, line));
    }
    const auto maximum = *std::max_element(crossings.begin(), crossings.end());
    const double average =
        static_cast<double>(std::accumulate(crossings.begin(), crossings.end(), std::size_t{})) /
        crossings.size();
    std::cout << "root-level crossings: max=" << maximum << ", avg=" << average
              << ", bound O(sqrt(r))=O(" << std::sqrt(64.0) << ")\n";
    return matches == 60 ? 0 : 1;
} catch (const std::exception &error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
}
