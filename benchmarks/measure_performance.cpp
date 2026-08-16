#include "matousek_partition_tree/core.hpp"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct Args {
    std::vector<std::size_t> sizes{200, 400, 800};
    std::size_t queries{100};
    int repeats{3};
    std::size_t r{25};
    std::size_t leaf_size{32};
    std::uint64_t seed{7};
    std::string json_path;
};

struct Row {
    std::size_t n{};
    std::size_t nodes{};
    double build_seconds{};
    double tree_us{};
    double brute_us{};
};

Args parse_args(int argc, char **argv) {
    Args args;
    for (int i = 1; i < argc; ++i) {
        const std::string option = argv[i];
        auto value = [&]() -> std::string {
            if (++i >= argc)
                throw std::invalid_argument("missing value for " + option);
            return argv[i];
        };
        if (option == "--smoke") {
            args.sizes = {120};
            args.queries = 8;
            args.repeats = 1;
            args.leaf_size = 16;
        } else if (option == "--sizes") {
            args.sizes.clear();
            while (i + 1 < argc && std::string(argv[i + 1]).rfind("--", 0) != 0) {
                args.sizes.push_back(std::stoull(argv[++i]));
            }
            if (args.sizes.empty())
                throw std::invalid_argument("--sizes requires values");
        } else if (option == "--queries") {
            args.queries = std::stoull(value());
        } else if (option == "--repeats") {
            args.repeats = std::stoi(value());
        } else if (option == "--r") {
            args.r = std::stoull(value());
        } else if (option == "--leaf-size") {
            args.leaf_size = std::stoull(value());
        } else if (option == "--seed") {
            args.seed = std::stoull(value());
        } else if (option == "--json") {
            args.json_path = value();
        } else {
            throw std::invalid_argument("unknown option: " + option);
        }
    }
    if (args.queries == 0 || args.repeats < 1 || args.r < 2 || args.leaf_size == 0) {
        throw std::invalid_argument("queries/repeats/leaf-size must be positive and r >= 2");
    }
    if (std::any_of(args.sizes.begin(), args.sizes.end(),
                    [&](auto n) { return n <= args.leaf_size; })) {
        throw std::invalid_argument("every size must exceed leaf-size");
    }
    return args;
}

std::size_t node_count(const mpt::PNode &node) {
    std::size_t count = 1;
    for (const auto &child : node.children)
        count += node_count(*child);
    return count;
}

template <class Function>
double median_microseconds(Function function, int repeats, std::size_t item_count) {
    std::vector<double> samples;
    for (int i = 0; i < repeats; ++i) {
        const auto start = std::chrono::steady_clock::now();
        function();
        samples.push_back(
            std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - start)
                .count() /
            item_count);
    }
    std::sort(samples.begin(), samples.end());
    return samples[samples.size() / 2];
}

Row measure(std::size_t n, const Args &args) {
    std::mt19937_64 rng(args.seed);
    std::uniform_int_distribution<std::int64_t> point_distribution(0, 10'000);
    std::uniform_int_distribution<std::int64_t> query_distribution(-10'000, 10'000);
    std::vector<mpt::Point> points;
    std::vector<mpt::Halfplane> queries;
    for (std::size_t i = 0; i < n; ++i) {
        points.push_back({mpt::rational(point_distribution(rng), 10'000),
                          mpt::rational(point_distribution(rng), 10'000)});
    }
    while (queries.size() < args.queries) {
        mpt::Halfplane query{mpt::rational(query_distribution(rng), 10'000),
                             mpt::rational(query_distribution(rng), 10'000),
                             mpt::rational(query_distribution(rng), 10'000)};
        if (query.a != 0 || query.b != 0)
            queries.push_back(std::move(query));
    }

    const auto build_start = std::chrono::steady_clock::now();
    auto tree = mpt::build_tree(points, args.r, args.leaf_size, args.seed);
    const double build_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - build_start).count();
    std::vector<std::size_t> expected;
    std::vector<std::size_t> actual;
    const double brute_us = median_microseconds(
        [&] {
            expected.clear();
            for (const auto &query : queries) {
                expected.push_back(std::count_if(points.begin(), points.end(), [&](const auto &p) {
                    return mpt::halfplane_side(query, p) >= 0;
                }));
            }
        },
        args.repeats, args.queries);
    const double tree_us = median_microseconds(
        [&] {
            actual.clear();
            for (const auto &query : queries)
                actual.push_back(mpt::query_count(*tree, query));
        },
        args.repeats, args.queries);
    if (actual != expected)
        throw std::logic_error("tree query disagrees with brute force");
    return {n, node_count(*tree), build_seconds, tree_us, brute_us};
}
} // namespace

int main(int argc, char **argv) try {
    const Args args = parse_args(argc, argv);
    std::vector<Row> rows;
    std::cout
        << "| n | nodes | build (s) | tree query (us) | brute query (us) | speedup | exact |\n"
           "|---:|---:|---:|---:|---:|---:|---:|\n";
    for (const auto n : args.sizes) {
        rows.push_back(measure(n, args));
        const auto &row = rows.back();
        std::cout << "| " << n << " | " << row.nodes << " | " << std::fixed << std::setprecision(3)
                  << row.build_seconds << " | " << std::setprecision(1) << row.tree_us << " | "
                  << row.brute_us << " | " << std::setprecision(2) << row.brute_us / row.tree_us
                  << "x | " << args.queries << '/' << args.queries << " |\n";
    }
    if (!args.json_path.empty()) {
        std::ofstream output(args.json_path);
        if (!output)
            throw std::runtime_error("cannot open JSON output: " + args.json_path);
        output << "{\n  \"seed\": " << args.seed << ",\n  \"repeats\": " << args.repeats
               << ",\n  \"results\": [\n";
        for (std::size_t i = 0; i < rows.size(); ++i) {
            const auto &row = rows[i];
            output << "    {\"n\": " << row.n << ", \"r\": " << args.r
                   << ", \"leaf_size\": " << args.leaf_size << ", \"queries\": " << args.queries
                   << ", \"tree_nodes\": " << row.nodes
                   << ", \"build_seconds\": " << row.build_seconds
                   << ", \"tree_query_us\": " << row.tree_us
                   << ", \"brute_query_us\": " << row.brute_us
                   << ", \"query_speedup\": " << row.brute_us / row.tree_us
                   << ", \"exact_matches\": " << args.queries << "}"
                   << (i + 1 == rows.size() ? "\n" : ",\n");
        }
        output << "  ]\n}\n";
    }
    return 0;
} catch (const std::exception &error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
}
