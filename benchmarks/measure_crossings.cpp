#include "matousek_partition_tree/core.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr std::array<std::size_t, 3> nominal_r_values{25, 36, 64};
constexpr std::size_t query_line_count = 200;

struct Row {
    std::size_t nominal_r{};
    std::size_t group_size{};
    double actual_ratio{};
    std::size_t groups{};
    std::size_t q_size{};
    int test_set_max{};
    std::size_t random_sample_max{};
};

struct BuiltPartition {
    std::size_t nominal_r{};
    std::size_t group_size{};
    mpt::PartitionResult result;
};

std::vector<mpt::Point> make_points(std::size_t n, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<std::int64_t> distribution(0, 10'000);
    std::vector<mpt::Point> points;
    points.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        points.push_back(
            {mpt::rational(distribution(rng), 10'000), mpt::rational(distribution(rng), 10'000)});
    }
    return points;
}

std::vector<mpt::Line> make_query_lines(std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<std::int64_t> distribution(-20'000, 20'000);
    std::vector<mpt::Line> lines;
    lines.reserve(query_line_count);
    for (std::size_t i = 0; i < query_line_count; ++i) {
        lines.push_back(
            {mpt::rational(distribution(rng), 10'000), mpt::rational(distribution(rng), 10'000)});
    }
    return lines;
}

std::pair<std::uint64_t, std::vector<BuiltPartition>>
build_with_common_seed(const std::vector<mpt::Point> &points, std::uint64_t initial_seed) {
    std::string last_error;
    for (std::uint64_t offset = 0; offset < 64; ++offset) {
        const std::uint64_t construction_seed = initial_seed + offset;
        std::vector<BuiltPartition> built;
        try {
            for (const auto nominal_r : nominal_r_values) {
                const std::size_t group_size = points.size() / nominal_r;
                std::mt19937_64 rng(construction_seed);
                built.push_back(
                    {nominal_r, group_size, mpt::simplicial_partition(points, group_size, rng)});
            }
            return {construction_seed, std::move(built)};
        } catch (const std::runtime_error &error) {
            last_error = error.what();
        }
    }
    throw std::runtime_error("no common construction seed succeeded: " + last_error);
}

void write_plot(const std::string &path, const std::vector<Row> &rows, std::size_t n,
                std::uint64_t point_seed, std::uint64_t construction_seed,
                std::uint64_t query_seed) {
    std::ofstream svg(path);
    if (!svg)
        throw std::runtime_error("cannot open plot output: " + path);

    constexpr double width = 820;
    constexpr double height = 580;
    constexpr double left = 92;
    constexpr double right = 30;
    constexpr double top = 112;
    constexpr double bottom = 82;
    constexpr double frame_width = width - left - right;
    constexpr double frame_height = height - top - bottom;
    const auto maximum = std::max_element(rows.begin(), rows.end(), [](const Row &a, const Row &b) {
        return std::max({a.groups, static_cast<std::size_t>(a.test_set_max), a.random_sample_max}) <
               std::max({b.groups, static_cast<std::size_t>(b.test_set_max), b.random_sample_max});
    });
    const double max_value =
        std::max({maximum->groups, static_cast<std::size_t>(maximum->test_set_max),
                  maximum->random_sample_max});
    const double y_max = std::ceil(max_value * 1.08 / 10.0) * 10.0;
    const auto x = [&](double value) {
        constexpr double inset = 18;
        return left + inset +
               (value - nominal_r_values.front()) /
                   static_cast<double>(nominal_r_values.back() - nominal_r_values.front()) *
                   (frame_width - 2 * inset);
    };
    const auto y = [&](double value) { return top + frame_height * (1.0 - value / y_max); };

    svg << std::fixed << std::setprecision(2);
    svg << "<svg xmlns='http://www.w3.org/2000/svg' width='" << width << "' height='" << height
        << "' viewBox='0 0 " << width << ' ' << height << "'>"
        << "<rect width='" << width << "' height='" << height << "' fill='#ffffff'/>"
        << "<style>text{font-family:system-ui,-apple-system,sans-serif;fill:#202124}"
           ".tick{font-size:12px}.label{font-size:12px;font-weight:500}"
           ".legend{font-size:13px}.meta{font-size:12px;fill:#5f6368}</style>"
        << "<text x='" << width / 2
        << "' y='30' text-anchor='middle' font-size='21' font-weight='500'>Measured crossing counts"
           "</text>"
        << "<text class='meta' x='" << width / 2 << "' y='51' text-anchor='middle'>n=" << n
        << ", point seed=" << point_seed << ", construction seed=" << construction_seed
        << ", query-line seed=" << query_seed << "</text>"
        << "<text class='meta' x='" << width / 2
        << "' y='69' text-anchor='middle'>Random-sample maximum is over the same "
        << query_line_count << " query lines for every partition.</text>";

    for (int tick = 0; tick <= static_cast<int>(y_max); tick += 10) {
        const double tick_y = y(tick);
        svg << "<line x1='" << left << "' y1='" << tick_y << "' x2='" << left + frame_width
            << "' y2='" << tick_y << "' stroke='" << (tick == 0 ? "#5f6368" : "#e0e0e0")
            << "' stroke-width='" << (tick == 0 ? 1.2 : 1.0) << "'/>"
            << "<text class='tick' x='" << left - 12 << "' y='" << tick_y + 4
            << "' text-anchor='end'>" << tick << "</text>";
    }
    svg << "<rect x='" << left << "' y='" << top << "' width='" << frame_width << "' height='"
        << frame_height << "' fill='none' stroke='#9aa0a6'/>";

    for (const auto &row : rows) {
        const double tick_x = x(row.nominal_r);
        svg << "<line x1='" << tick_x << "' y1='" << top + frame_height << "' x2='" << tick_x
            << "' y2='" << top + frame_height + 6 << "' stroke='#5f6368'/>"
            << "<text class='tick' x='" << tick_x << "' y='" << top + frame_height + 23
            << "' text-anchor='middle'>" << row.nominal_r << "</text>";
    }

    struct Series {
        const char *name;
        const char *color;
        std::vector<double> values;
        int label_offset;
        bool square_marker;
        const char *dash;
    };
    const std::array<Series, 3> series{{
        {"Test-set maximum", "#3f6fa9", {}, 20, false, ""},
        {"Random-sample maximum", "#d64f4f", {}, 20, false, ""},
        {"Group count", "#70757a", {}, -12, true, "6 4"},
    }};
    std::array<Series, 3> plotted = series;
    for (const auto &row : rows) {
        plotted[0].values.push_back(row.test_set_max);
        plotted[1].values.push_back(row.random_sample_max);
        plotted[2].values.push_back(row.groups);
    }
    for (const auto &item : plotted) {
        svg << "<polyline fill='none' stroke='" << item.color
            << "' stroke-width='2.5' stroke-linejoin='round' stroke-dasharray='" << item.dash
            << "' points='";
        for (std::size_t i = 0; i < rows.size(); ++i)
            svg << x(rows[i].nominal_r) << ',' << y(item.values[i]) << ' ';
        svg << "'/>";
        for (std::size_t i = 0; i < rows.size(); ++i) {
            const double point_x = x(rows[i].nominal_r);
            const double point_y = y(item.values[i]);
            if (item.square_marker) {
                svg << "<rect x='" << point_x - 5 << "' y='" << point_y - 5
                    << "' width='10' height='10' fill='none' stroke='" << item.color
                    << "' stroke-width='2.5'/>";
            } else {
                svg << "<circle cx='" << point_x << "' cy='" << point_y
                    << "' r='4.5' fill='#ffffff' stroke='" << item.color
                    << "' stroke-width='2.5'/>";
            }
            svg << "<text class='label' x='" << point_x << "' y='" << point_y + item.label_offset
                << "' text-anchor='middle' style='fill:" << item.color << "'>"
                << static_cast<std::size_t>(item.values[i]) << "</text>";
        }
    }

    constexpr std::array<double, 3> legend_x{170, 365, 595};
    for (std::size_t i = 0; i < plotted.size(); ++i) {
        svg << "<line x1='" << legend_x[i] << "' y1='91' x2='" << legend_x[i] + 24
            << "' y2='91' stroke='" << plotted[i].color << "' stroke-width='2.5' "
            << "stroke-dasharray='" << plotted[i].dash << "'/>";
        if (plotted[i].square_marker) {
            svg << "<rect x='" << legend_x[i] + 8 << "' y='87' width='8' height='8' fill='none' "
                << "stroke='" << plotted[i].color << "' stroke-width='2'/>";
        } else {
            svg << "<circle cx='" << legend_x[i] + 12 << "' cy='91' r='3.5' fill='#fff' stroke='"
                << plotted[i].color << "' stroke-width='2'/>";
        }
        svg << "<text class='legend' x='" << legend_x[i] + 31 << "' y='95'>" << plotted[i].name
            << "</text>";
    }

    svg << "<text x='" << left + frame_width / 2 << "' y='" << height - 32
        << "' text-anchor='middle' font-size='13'>Nominal partition parameter R</text>"
        << "<text x='22' y='" << top + frame_height / 2
        << "' text-anchor='middle' font-size='13' transform='rotate(-90 22 "
        << top + frame_height / 2 << ")'>Count</text>"
        << "<text class='meta' x='" << left + frame_width / 2 << "' y='" << height - 12
        << "' text-anchor='middle'>Group size uses integer division; group count can differ from "
           "nominal R.</text>"
        << "</svg>\n";
}
} // namespace

int main(int argc, char **argv) try {
    std::vector<std::string> positional;
    std::string plot_path;
    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        if (argument == "--plot") {
            if (i + 1 >= argc)
                throw std::invalid_argument("--plot requires an output path");
            plot_path = argv[++i];
        } else if (argument.rfind("--", 0) == 0) {
            throw std::invalid_argument("unknown option: " + argument);
        } else {
            positional.push_back(argument);
        }
    }
    if (positional.size() > 2)
        throw std::invalid_argument("usage: measure-crossings [n] [seed] [--plot out.svg]");

    const std::size_t n = positional.empty() ? 1200 : std::stoull(positional[0]);
    const std::uint64_t point_seed = positional.size() < 2 ? 7 : std::stoull(positional[1]);
    if (n < nominal_r_values.back())
        throw std::invalid_argument("n must be at least 64");
    const std::uint64_t query_seed = point_seed + 1;
    const auto points = make_points(n, point_seed);
    const auto query_lines = make_query_lines(query_seed);
    auto [construction_seed, built] = build_with_common_seed(points, point_seed);

    std::vector<Row> rows;
    for (const auto &partition : built) {
        std::size_t random_sample_max = 0;
        for (const auto &line : query_lines) {
            const auto crossings = std::count_if(
                partition.result.groups.begin(), partition.result.groups.end(),
                [&](const auto &group) { return mpt::line_crosses_tri(line, group.simplex); });
            random_sample_max = std::max(random_sample_max, static_cast<std::size_t>(crossings));
        }
        rows.push_back({partition.nominal_r, partition.group_size,
                        static_cast<double>(n) / partition.group_size,
                        partition.result.groups.size(), partition.result.stats.Q.size(),
                        partition.result.stats.K_Q, random_sample_max});
    }

    std::cout << "n=" << n << ", point seed=" << point_seed
              << ", construction seed=" << construction_seed << ", query-line seed=" << query_seed
              << ", shared query lines=" << query_line_count << "\n\n"
              << "| nominal R | group size | actual n/s | groups | |Q| | test-set max | "
                 "random-sample max |\n"
              << "|---:|---:|---:|---:|---:|---:|---:|\n";
    for (const auto &row : rows) {
        std::cout << "| " << row.nominal_r << " | " << row.group_size << " | " << std::fixed
                  << std::setprecision(2) << row.actual_ratio << " | " << row.groups << " | "
                  << row.q_size << " | " << row.test_set_max << " | " << row.random_sample_max
                  << " |\n";
    }
    if (!plot_path.empty()) {
        write_plot(plot_path, rows, n, point_seed, construction_seed, query_seed);
        std::cout << "\nplot written: " << plot_path << '\n';
    }
    return 0;
} catch (const std::exception &error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
}
