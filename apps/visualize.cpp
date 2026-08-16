#include "matousek_partition_tree/core.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <optional>
#include <random>
#include <string>
#include <vector>

namespace {
double px(const mpt::Rational &x, double offset) {
    return offset + 35.0 + 430.0 * mpt::to_double(x);
}
double py(const mpt::Rational &y) { return 485.0 - 430.0 * mpt::to_double(y); }
} // namespace

int main(int argc, char **argv) try {
    const std::size_t n = argc > 1 ? std::stoull(argv[1]) : 1200;
    const std::uint64_t seed = argc > 2 ? std::stoull(argv[2]) : 42;
    const std::string path = argc > 3 ? argv[3] : "assets/partition_tree_example.svg";
    const std::size_t nominal_r = argc > 4 ? std::stoull(argv[4]) : 64;
    const double cutting_scale = argc > 5 ? std::stod(argv[5]) : 0.35;
    if (nominal_r < 2 || n < nominal_r)
        throw std::invalid_argument("require n >= nominal R >= 2");
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<std::int64_t> distribution(0, 10'000);
    std::vector<mpt::Point> points;
    for (std::size_t i = 0; i < n; ++i) {
        points.push_back(
            {mpt::rational(distribution(rng), 10'000), mpt::rational(distribution(rng), 10'000)});
    }
    const std::size_t group_size = n / nominal_r;
    std::optional<mpt::PartitionResult> built;
    std::uint64_t construction_seed = seed;
    std::string last_error;
    for (std::uint64_t offset = 0; offset < 64 && !built; ++offset) {
        construction_seed = seed + offset;
        try {
            std::mt19937_64 build_rng(construction_seed);
            built = mpt::simplicial_partition(points, group_size, build_rng, cutting_scale);
        } catch (const std::runtime_error &error) {
            last_error = error.what();
        }
    }
    if (!built)
        throw std::runtime_error("all deterministic construction attempts failed: " + last_error);
    auto &result = *built;
    const mpt::Line query{mpt::rational(3, 10), mpt::rational(35, 100)};
    std::ofstream svg(path);
    if (!svg)
        throw std::runtime_error("cannot open output file: " + path);
    svg << "<svg xmlns='http://www.w3.org/2000/svg' width='1500' height='570' viewBox='0 0 1500 "
           "570'>"
           "<rect x='0' y='0' width='1500' height='570' fill='#ffffff'/>"
           "<defs><clipPath id='p0'><rect x='35' y='55' width='430' height='430'/></clipPath>"
           "<clipPath id='p1'><rect x='535' y='55' width='430' height='430'/></clipPath>"
           "<clipPath id='p2'><rect x='1035' y='55' width='430' height='430'/></clipPath></defs>"
           "<style>text{font-family:system-ui,sans-serif;fill:#111}.title{font-size:17px;font-"
           "weight:600}</style>";
    const char *colors[] = {"#4e79a7", "#f28e2b", "#e15759", "#76b7b2",
                            "#59a14f", "#edc949", "#af7aa1", "#ff9da7"};
    const double background_radius = n >= 5'000 ? 0.55 : 1.2;
    const double selected_radius = n >= 5'000 ? 1.8 : 3.0;
    const double crossed_radius = n >= 5'000 ? 1.0 : 2.6;
    std::vector<std::size_t> selected;
    for (std::size_t i = 0; i < 8; ++i)
        selected.push_back(i * (result.groups.size() - 1) / 7);
    int panel = 0;
    for (const double offset : {0.0, 500.0, 1000.0}) {
        svg << "<rect x='" << offset + 35
            << "' y='55' width='430' height='430' fill='none' stroke='#bbb'/>";
        for (const auto &point : points)
            svg << "<circle clip-path='url(#p" << panel << ")' cx='" << px(point.x, offset)
                << "' cy='" << py(point.y) << "' r='" << background_radius << "' fill='#ccc'/>";
        ++panel;
    }
    for (std::size_t j = 0; j < selected.size(); ++j) {
        const auto &group = result.groups[selected[j]];
        for (const auto &point : group.points) {
            svg << "<circle clip-path='url(#p0)' cx='" << px(point.x, 0) << "' cy='" << py(point.y)
                << "' r='" << selected_radius << "' fill='" << colors[j] << "'/>";
            svg << "<circle clip-path='url(#p1)' cx='" << px(point.x, 500) << "' cy='"
                << py(point.y) << "' r='" << selected_radius << "' fill='" << colors[j] << "'/>";
        }
        svg << "<polygon clip-path='url(#p1)' points='";
        for (const auto &vertex : group.simplex)
            svg << px(vertex.x, 500) << ',' << py(vertex.y) << ' ';
        svg << "' fill='none' stroke='" << colors[j] << "' stroke-width='1.5'/>";
    }
    std::size_t crossed = 0;
    for (const auto &group : result.groups) {
        const bool crosses = mpt::line_crosses_tri(query, group.simplex);
        crossed += crosses;
        for (const auto &point : group.points) {
            svg << "<circle clip-path='url(#p2)' cx='" << px(point.x, 1000) << "' cy='"
                << py(point.y) << "' r='" << (crosses ? crossed_radius : background_radius)
                << "' fill='" << (crosses ? "#d85a30" : "#bbb") << "'/>";
        }
    }
    svg << "<line clip-path='url(#p2)' x1='1035' y1='" << py(mpt::rational(35, 100))
        << "' x2='1465' y2='" << py(mpt::rational(65, 100))
        << "' stroke='black' stroke-width='2'/>";
    svg << "<text class='title' x='250' y='30' text-anchor='middle'>Selected point groups</text>"
           "<text class='title' x='750' y='30' text-anchor='middle'>Groups and containing "
           "simplices</text>"
           "<text class='title' x='1250' y='30' text-anchor='middle'>Query-line crossings</text>"
           "<text x='750' y='548' text-anchor='middle' font-size='13'>";
    svg << "n=" << n << ", nominal R=" << nominal_r << ", group size=" << group_size
        << ", groups=" << result.groups.size() << ", construction seed=" << construction_seed
        << ", cutting scale=" << cutting_scale << ", test-set maximum=" << result.stats.K_Q
        << ", query crossings=" << crossed << '/' << result.groups.size()
        << "</text>"
           "</svg>\n";
    std::cout << "groups=" << result.groups.size() << ", K_Q=" << result.stats.K_Q
              << ", crossed=" << crossed << ", written: " << path << '\n';
    return 0;
} catch (const std::exception &error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
}
