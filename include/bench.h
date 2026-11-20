// bench.h - every benchmark behind the report's complexity, comparison and improvement sections
//
// Timing protocol: 1 warm-up run, then 5 timed runs; the CSVs report the median and the minimum.
// Per-seam timings take the energy map as given (computed outside the timer). Built with -O2.
// Each CSV starts with a "# machine" comment line. 
//
// Suites (CSV in --out):
//   scaling  bench_scaling.csv  per-seam time vs n (square random energy) for every finder
//   carve    bench_carve.csv    full 50% carve of each images/* at input width 128 and 320
//   flow     bench_flow.csv     FF vs EK vs Dinic on the seam graph
//   improve  bench_improve.csv  memory, energy update modes, threads, max-flow solvers
//   sobel    bench_sobel.csv    direct vs separable vs D&C vs threaded energy map (Q2b)
//   count    seam_counts.csv    number of seams vs 2^(h-1) and w 3^(h-1) (Q2a)
#pragma once

#include <algorithm>
#include <thread>

#include "carve.h"

namespace sc {

struct BenchArgs {
    std::string suite = "all", out = "results";
    bool quick = false;
};

namespace bench {

constexpr double kBudgetMs = 10000;

inline std::string machine_note() {
    std::string cpu;
    if (const char* id = std::getenv("PROCESSOR_IDENTIFIER")) cpu = id;  // Windows
    std::FILE* f = std::fopen("/proc/cpuinfo", "r");                     // Linux
    if (f) {
        char line[256];
        while (cpu.empty() && std::fgets(line, sizeof line, f)) {
            std::string l = line;
            if (l.rfind("model name", 0) == 0) cpu = l.substr(l.find(':') + 2), cpu.pop_back();
        }
        std::fclose(f);
    }
    cpu.erase(std::remove(cpu.begin(), cpu.end(), ','), cpu.end());  // the note itself is comma separated
    return "cpu=" + (cpu.empty() ? std::string("unknown") : cpu) +
           ", threads=" + std::to_string(std::thread::hardware_concurrency()) + ", -O2";
}

struct Timed {
    bool skipped = false;
    double median = 0, min = 0;
};

// 1 warm-up + `runs` timed runs. Skipped (no timed runs) when the warm-up exceeds the budget.
inline Timed timed(const std::function<void()>& f, int runs = 5) {
    Timer t;
    f();
    Timed r;
    if (t.ms() > kBudgetMs) {
        r.skipped = true;
        return r;
    }
    std::vector<double> v;
    for (int i = 0; i < runs; ++i) {
        t.start();
        f();
        v.push_back(t.ms());
    }
    r.median = median(v);
    r.min = *std::min_element(v.begin(), v.end());
    return r;
}

inline void progress(const std::string& suite, const std::string& what, const Timed& t) {
    if (t.skipped)
        std::printf("[%s] %s skipped (warm-up over %.0f s)\n", suite.c_str(), what.c_str(), kBudgetMs / 1000);
    else
        std::printf("[%s] %s %.3f ms\n", suite.c_str(), what.c_str(), t.median);
    std::fflush(stdout);
}

inline std::vector<std::string> image_files() {
    std::vector<std::string> files;
    if (!std::filesystem::exists("images")) return files;
    for (const auto& entry : std::filesystem::directory_iterator("images")) {
        std::string ext = entry.path().extension().string();
        for (char& ch : ext) ch = char(std::tolower(static_cast<unsigned char>(ch)));
        if (ext == ".png" || ext == ".jpg" || ext == ".jpeg") files.push_back(entry.path().generic_string());
    }
    std::sort(files.begin(), files.end());
    return files;
}

inline Image resize_to_width(const Image& img, int w) {
    if (img.w <= w) return img;
    return resize_nearest(img, w, std::max(1, int(std::lround(double(img.h) * w / img.w))));
}

// Total energy actually removed, measured on the true (freshly computed) energy at each step.
// The carve's own sum is not comparable for UpdateMode::Stale, whose map is out of date.
inline int64_t true_removed_energy(const Image& in, const std::vector<std::vector<int32_t>>& removed, EnergyKind k) {
    Image img = in;
    IndexGrid idx(in.h, in.w);
    for (int i = 0; i < int(idx.a.size()); ++i) idx.a[i] = i;
    int64_t total = 0;
    for (const auto& ids : removed) {
        Energy e = compute_energy(to_gray(img), k);
        Seam s(img.h);
        for (int r = 0; r < img.h; ++r) {
            const int32_t* row = &idx.at(r, 0);
            s[r] = int(std::find(row, row + img.w, ids[r]) - row);
        }
        total += seam_cost(e, s);
        img = remove_seam(img, s);
        idx = remove_seam(idx, s);
    }
    return total;
}

inline std::string fmt_ms(const Timed& t) {
    if (t.skipped) return "skipped";
    char buf[32];
    std::snprintf(buf, sizeof buf, "%.6g", t.median);
    return buf;
}
inline std::string fmt_min(const Timed& t) {
    if (t.skipped) return "";
    char buf[32];
    std::snprintf(buf, sizeof buf, "%.6g", t.min);
    return buf;
}

inline void suite_scaling(const BenchArgs& a) {
    Csv csv(a.out + "/bench_scaling.csv", {"method", "n", "h", "w", "median_ms", "min_ms", "cells", "augment", "bfs_visits", "phases"}, machine_note());
    const std::vector<int> dp_sizes = {64, 128, 256, 512, 1024}, memo_sizes = {64, 128, 256, 512},
                           flow_sizes = {16, 32, 48, 64, 96, 128};
    const std::vector<std::pair<std::string, std::vector<int>>> specs = {
        {"dp", dp_sizes},  {"rolling", dp_sizes}, {"memo", memo_sizes}, {"mt", dp_sizes},
        {"greedy", dp_sizes}, {"ff", flow_sizes}, {"ek", flow_sizes},   {"dinic", flow_sizes}};
    for (const auto& [method, all_sizes] : specs) {
        std::vector<int> sizes(all_sizes.begin(), all_sizes.end() - (a.quick ? 2 : 0));
        Finder f = get_finder(method, 4);
        bool over_budget = false;
        for (int n : sizes) {
            if (over_budget) {
                csv.row(method, n, n, n, std::string("skipped"), std::string(""), 0, 0, 0, 0);
                continue;
            }
            std::mt19937 rng(1000 + n);
            Energy e = random_energy(rng, n, n, 2040);
            Stats st;
            Timed t = timed([&] {
                st = Stats{};
                f(e, st);
            });
            over_budget = t.skipped;
            progress("scaling", method + " n=" + std::to_string(n), t);
            csv.row(method, n, n, n, fmt_ms(t), fmt_min(t), st.cells, st.augment, st.bfs_visits, st.phases);
        }
    }
}

inline void suite_carve(const BenchArgs& a) {
    std::vector<std::string> files = image_files();
    if (files.empty()) {
        std::printf("[carve] warning: no images in images/, suite skipped\n");
        return;
    }
    Csv csv(a.out + "/bench_carve.csv", {"image", "method", "w0", "h0", "W", "H", "total_ms", "removed_energy", "mean_seam_cost"}, machine_note() + ", single run per row");
    std::vector<int> widths = a.quick ? std::vector<int>{128} : std::vector<int>{128, 320};
    for (const std::string& file : files) {
        std::string stem = std::filesystem::path(file).stem().string();
        Image full = load_image(file);
        for (int w0 : widths) {
            Image in = resize_to_width(full, w0);
            if (in.w != w0) continue;  // image narrower than this input width
            const int W = w0 / 2, H = in.h;
            std::vector<std::string> methods = {"dp", "greedy"};
            if (w0 <= 128) methods.push_back("dinic");  // graph cut only on the small input: it takes minutes at 320 px
            for (const std::string& m : methods) {
                CarveOptions o;
                o.method = m;
                CarveResult r = carve(in, W, H, o);
                double mean = r.seam_costs.empty() ? 0 : double(r.removed_energy) / r.seam_costs.size();
                csv.row(stem, m, w0, in.h, W, H, r.ms, r.removed_energy, mean);
                std::printf("[carve] %s w0=%d %s %.1f ms, removed energy %lld\n", stem.c_str(), w0, m.c_str(), r.ms,
                            (long long)r.removed_energy);
                std::fflush(stdout);
            }
            for (const std::string& k : baseline_names()) {
                Timer t;
                Image out = run_baseline(k, in, W, H, EnergyKind::L1);
                double ms = t.ms();
                csv.row(stem, k, w0, in.h, W, H, ms, std::string(""), std::string(""));
            }
        }
    }
}

inline void suite_flow(const BenchArgs& a) {
    Csv csv(a.out + "/bench_flow.csv", {"solver", "n", "flow", "median_ms", "augment", "bfs_visits", "phases"}, machine_note() + ", each timed run copies the prebuilt graph (O(m), small next to the solve)");
    std::vector<int> sizes = {8, 16, 24, 32, 48, 64};
    if (a.quick) sizes.resize(4);
    using MaxFlow = int64_t (*)(FlowGraph&, Stats&);
    const std::pair<const char*, MaxFlow> solvers[] = {{"ff", maxflow_ff}, {"ek", maxflow_ek}, {"dinic", maxflow_dinic}};
    for (int n : sizes) {
        std::mt19937 rng(2000 + n);
        SeamGraph base = build_seam_graph(random_energy(rng, n, n, 2040));
        for (const auto& [name, f] : solvers) {
            const MaxFlow solve = f;  // a lambda may not capture a structured binding in C++17
            Stats st;
            int64_t flow = 0;
            Timed t = timed([&] {
                FlowGraph g = base.g;
                st = Stats{};
                flow = solve(g, st);
            });
            progress("flow", std::string(name) + " n=" + std::to_string(n), t);
            csv.row(std::string(name), n, flow, fmt_ms(t), st.augment, st.bfs_visits, st.phases);
        }
    }
}

inline void suite_improve(const BenchArgs& a) {
    Csv csv(a.out + "/bench_improve.csv", {"experiment", "variant", "param", "median_ms", "bytes", "removed_energy"}, machine_note());
    const std::string none;

    // (a) Do we need the full DP table? Time and memory of full vs rolling rows.
    {
        const int n = a.quick ? 512 : 1024;
        std::mt19937 rng(3001);
        Energy e = random_energy(rng, n, n, 2040);
        Stats st;
        const long long full_bytes = 8LL * n * (n + 2) + 1LL * n * n, rolling_bytes = 2 * 8LL * (n + 2) + 1LL * n * n;
        Timed full = timed([&] { dp_seam(e, st); }), roll = timed([&] { dp_seam_rolling(e, st); });
        progress("improve", "full table n=" + std::to_string(n), full);
        progress("improve", "rolling n=" + std::to_string(n), roll);
        csv.row(std::string("a_dp_memory"), std::string("full table"), "n=" + std::to_string(n), fmt_ms(full), full_bytes, none);
        csv.row(std::string("a_dp_memory"), std::string("rolling rows"), "n=" + std::to_string(n), fmt_ms(roll), rolling_bytes, none);
    }
    // (b) Energy update after each removal: full recompute vs local band vs stale (control).
    {
        std::vector<std::string> files = image_files();
        const int w0 = a.quick ? 160 : 320;
        Image in;
        std::string label;
        if (!files.empty()) {
            in = resize_to_width(load_image(files.front()), w0);
            label = std::filesystem::path(files.front()).stem().string();
        } else {
            std::mt19937 rng(3002);
            in = random_image(rng, w0 * 3 / 4, w0);
            label = "random";
        }
        for (UpdateMode m : {UpdateMode::Full, UpdateMode::Local, UpdateMode::Stale}) {
            CarveOptions o;
            o.update = m;
            CarveResult r;
            Timed t = timed([&] { r = carve(in, in.w / 2, in.h, o); }, 3);
            int64_t removed = true_removed_energy(in, r.removed_v, o.energy);
            progress("improve", "update " + to_string(m), t);
            csv.row(std::string("b_energy_update"), to_string(m), label + " w=" + std::to_string(in.w), fmt_ms(t), none, removed);
        }
    }
    // (c) Threads: the energy map parallelises well; the row-by-row DP does not.
    {
        const int ne = a.quick ? 1024 : 2048, nd = a.quick ? 512 : 1024;
        std::mt19937 rng(3003);
        Gray g(ne, ne);
        for (auto& v : g.a) v = uint8_t(rng());
        Energy e = random_energy(rng, nd, nd, 2040);
        for (int t : {1, 2, 4, 8}) {
            Timed te = timed([&] { compute_energy(g, EnergyKind::L1, t); });
            progress("improve", "energy threads=" + std::to_string(t), te);
            csv.row(std::string("c_threads_energy"), "threads=" + std::to_string(t), "n=" + std::to_string(ne), fmt_ms(te), none, none);
        }
        for (int t : {1, 2, 4, 8}) {
            Stats st;
            Timed td = timed([&] {
                if (t == 1)
                    dp_seam(e, st);
                else
                    dp_seam_mt(e, st, t);
            });
            progress("improve", "dp threads=" + std::to_string(t), td);
            csv.row(std::string("c_threads_dp"), t == 1 ? std::string("serial dp") : "threads=" + std::to_string(t),
                    "n=" + std::to_string(nd), fmt_ms(td), none, none);
        }
    }
    // (d) Better max-flow algorithms for the graph-cut seam.
    {
        const int n = a.quick ? 48 : 64;
        std::mt19937 rng(3004);
        Energy e = random_energy(rng, n, n, 2040);
        for (std::string m : {"ff", "ek", "dinic"}) {
            Finder f = get_finder(m);
            Stats st;
            Timed t = timed([&] { f(e, st); });
            progress("improve", "maxflow " + m, t);
            csv.row(std::string("d_maxflow"), m, "n=" + std::to_string(n), fmt_ms(t), none, none);
        }
    }
}

inline void suite_sobel(const BenchArgs& a) {
    Csv csv(a.out + "/bench_sobel.csv", {"variant", "n", "median_ms"}, machine_note());
    std::vector<int> sizes = {256, 512, 1024, 2048};
    if (a.quick) sizes.resize(2);
    for (int n : sizes) {
        std::mt19937 rng(4000 + n);
        Gray g(n, n);
        for (auto& v : g.a) v = uint8_t(rng());
        const std::pair<const char*, std::function<void()>> variants[] = {
            {"direct", [&] { sobel_direct(g, EnergyKind::L1); }},
            {"separable", [&] { sobel_separable(g, EnergyKind::L1); }},
            {"dc", [&] { sobel_dc(g, EnergyKind::L1); }},
            {"separable 4 threads", [&] { compute_energy(g, EnergyKind::L1, 4); }}};
        for (const auto& [name, f] : variants) {
            Timed t = timed(f);
            progress("sobel", std::string(name) + " n=" + std::to_string(n), t);
            csv.row(std::string(name), n, fmt_ms(t));
        }
    }
}

inline void suite_count(const BenchArgs& a) {
    Csv csv(a.out + "/seam_counts.csv", {"h", "w", "count", "lower_2pow", "upper_3pow"});
    for (int w : {2, 3, 10, 100})
        for (int h = 1; h <= 30; ++h) csv.row(h, w, count_seams(h, w), uint64_t(1) << (h - 1), double(w) * std::pow(3.0, h - 1));
    std::printf("[count] seam_counts.csv written\n");
}

}  // namespace bench

inline int run_bench(const BenchArgs& a) {
    const std::vector<std::pair<std::string, void (*)(const BenchArgs&)>> suites = {
        {"count", bench::suite_count},   {"sobel", bench::suite_sobel}, {"flow", bench::suite_flow},
        {"scaling", bench::suite_scaling}, {"improve", bench::suite_improve}, {"carve", bench::suite_carve}};
    bool found = false;
    Timer total;
    for (const auto& [name, run] : suites)
        if (a.suite == "all" || a.suite == name) {
            found = true;
            run(a);
        }
    if (!found) throw std::invalid_argument("unknown suite '" + a.suite + "' (all, scaling, carve, flow, improve, sobel, count)");
    std::printf("bench done in %.1f s, CSVs in %s/\n", total.ms() / 1000, a.out.c_str());
    return 0;
}

}  // namespace sc
