// main.cpp - the `seamcarve` command-line tool. Owner: Vera.
// Every command is a thin layer over the modules; run `seamcarve --help` for the usage text.
// Exit codes: 0 success, 1 runtime error (bad file, invalid size, ...), 2 bad usage.
#include <map>
#include <mutex>
#include <set>

#include "bench.h"

using namespace sc;

namespace {

const char* kUsage = R"(seamcarve - content-aware image resizing by seam carving

usage:
  seamcarve carve IN [-W w] [-H h] [-m METHOD] [-e l1|l2|grad] [-u full|local|stale] [-t THREADS]
                     [-o OUT] [--seams OVERLAY] [--frames DIR [--every K]] [--csv FILE]
  seamcarve baseline IN [-W w] [-H h] --kind crop|scale|column|rowmin [-e l1|l2|grad] [-o OUT]
  seamcarve energy IN [-e l1|l2|grad] [-t THREADS] [-o OUT]
  seamcarve dumpdp (--demo | IN)
  seamcarve order IN --min-width K [-o seams.json]
  seamcarve bench [--suite all|scaling|carve|flow|improve|sobel|count] [--quick] [--out DIR]
  seamcarve demo [--out DIR]
  seamcarve methods
  seamcarve --help

sizes (-W, -H, --min-width) are pixels (300) or a percentage of the input (50%); the default is 100%.
methods: dp (default), rolling, memo, mt, greedy, brute, ff, ek, dinic.
defaults: -e l1, -u local, -t 1, -o carved.png (energy.png, baseline.png, seams.json).

examples:
  seamcarve carve images/starry_night.jpg -W 50% -o out.png --seams seams.png
  seamcarve carve images/synthetic.png -W 60% -m dinic -o dinic.png
  seamcarve baseline images/synthetic.png -W 50% --kind scale -o scaled.png
  seamcarve dumpdp --demo
  seamcarve bench --quick
)";

struct UsageError : std::runtime_error {
    using std::runtime_error::runtime_error;
};

// Positional arguments plus "-x VALUE" / "--name VALUE" options and "--flag" switches.
struct Args {
    std::vector<std::string> pos;
    std::map<std::string, std::string> opt;
    std::set<std::string> flags;

    Args(int argc, char** argv, const std::set<std::string>& with_value, const std::set<std::string>& switches) {
        for (int i = 2; i < argc; ++i) {
            std::string a = argv[i];
            if (with_value.count(a)) {
                if (i + 1 >= argc) throw UsageError(a + " needs a value");
                opt[a] = argv[++i];
            } else if (switches.count(a)) {
                flags.insert(a);
            } else if (a.size() > 1 && a[0] == '-' && !std::isdigit(static_cast<unsigned char>(a[1]))) {
                throw UsageError("unknown option " + a);
            } else {
                pos.push_back(a);
            }
        }
    }
    std::string get(const std::string& k, const std::string& def) const {
        auto it = opt.find(k);
        return it == opt.end() ? def : it->second;
    }
    std::string input(size_t n = 1) const {
        if (pos.size() != n) throw UsageError(n == 1 ? "expected exactly one input image" : "wrong number of arguments");
        return pos[0];
    }
};

// Runs a parser of an argument value (size, method, energy, ...) and reports a bad value as bad usage.
template <class F>
auto arg(F parse) -> decltype(parse()) {
    try {
        return parse();
    } catch (const std::invalid_argument& e) {
        throw UsageError(e.what());
    }
}

int parse_int(const std::string& s, const std::string& what, int lo) {
    size_t used = 0;
    int v = 0;
    try {
        v = std::stoi(s, &used);
    } catch (const std::exception&) {
        used = 0;
    }
    if (used != s.size() || v < lo) throw UsageError(what + " must be an integer >= " + std::to_string(lo) + ", got '" + s + "'");
    return v;
}

int cmd_carve(int argc, char** argv) {
    Args a(argc, argv, {"-W", "-H", "-m", "-e", "-u", "-t", "-o", "--seams", "--frames", "--every", "--csv"}, {});
    Image in = load_image(a.input());
    int W = arg([&] { return parse_size(a.get("-W", "100%"), in.w); }), H = arg([&] { return parse_size(a.get("-H", "100%"), in.h); });
    CarveOptions o;
    o.method = a.get("-m", "dp");
    arg([&] { return get_finder(o.method); });  // reject an unknown method before any work
    o.energy = arg([&] { return parse_energy(a.get("-e", "l1")); });
    o.update = arg([&] { return parse_update(a.get("-u", "local")); });
    o.threads = parse_int(a.get("-t", "1"), "-t", 1);
    o.frames_dir = a.get("--frames", "");
    o.frame_every = o.frames_dir.empty() ? 0 : parse_int(a.get("--every", "1"), "--every", 1);
    if (o.method == "brute" && double(in.w) * std::pow(3.0, in.h - 1) > 5e7)
        throw std::invalid_argument("method brute is only for tiny images (w * 3^(h-1) <= 5e7); try -m dp");

    CarveResult r = carve(in, W, H, o);
    std::string out = a.get("-o", "carved.png");
    save_png(out, r.out);
    if (a.opt.count("--seams")) save_png(a.opt["--seams"], overlay_seams(in, r));
    if (a.opt.count("--csv")) {
        Csv csv(a.opt["--csv"], {"step", "direction", "seam_cost"});
        int step = 0;
        for (size_t i = 0; i < r.seam_costs.size(); ++i)
            csv.row(++step, std::string(i < r.removed_v.size() ? "vertical" : "horizontal"), r.seam_costs[i]);
    }
    std::printf("%dx%d -> %dx%d: %zu vertical + %zu horizontal seams, removed energy %lld, %.1f ms\n", in.w, in.h, r.out.w,
                r.out.h, r.removed_v.size(), r.removed_h.size(), (long long)r.removed_energy, r.ms);
    std::printf("method %s, energy %s, update %s -> %s\n", o.method.c_str(), to_string(o.energy).c_str(),
                to_string(o.update).c_str(), out.c_str());
    return 0;
}

int cmd_baseline(int argc, char** argv) {
    Args a(argc, argv, {"-W", "-H", "--kind", "-e", "-o"}, {});
    Image in = load_image(a.input());
    const std::vector<std::string> kinds = baseline_names();
    if (std::find(kinds.begin(), kinds.end(), a.get("--kind", "")) == kinds.end())
        throw UsageError("baseline needs --kind crop|scale|column|rowmin");
    int W = arg([&] { return parse_size(a.get("-W", "100%"), in.w); }), H = arg([&] { return parse_size(a.get("-H", "100%"), in.h); });
    EnergyKind k = arg([&] { return parse_energy(a.get("-e", "l1")); });
    Timer t;
    Image out = run_baseline(a.opt["--kind"], in, W, H, k);
    std::string path = a.get("-o", "baseline.png");
    save_png(path, out);
    std::printf("%dx%d -> %dx%d (%s), %.1f ms -> %s\n", in.w, in.h, out.w, out.h, a.opt["--kind"].c_str(), t.ms(), path.c_str());
    return 0;
}

int cmd_energy(int argc, char** argv) {
    Args a(argc, argv, {"-e", "-t", "-o"}, {});
    Image in = load_image(a.input());
    EnergyKind k = arg([&] { return parse_energy(a.get("-e", "l1")); });
    Energy e = compute_energy(to_gray(in), k, parse_int(a.get("-t", "1"), "-t", 1));
    std::string path = a.get("-o", "energy.png");
    save_png(path, energy_view(e));
    std::printf("%dx%d energy map (max %d) -> %s\n", e.w, e.h, *std::max_element(e.a.begin(), e.a.end()), path.c_str());
    return 0;
}

int cmd_dumpdp(int argc, char** argv) {
    Args a(argc, argv, {}, {"--demo"});
    Energy e;
    if (a.flags.count("--demo")) {
        if (!a.pos.empty()) throw UsageError("dumpdp takes --demo or one image, not both");
        e = dp_demo_energy();
        std::printf("DP table demo: the hand-checked 4x5 example (test DP6)\n\n");
    } else {
        Image in = load_image(a.input());
        if (in.w > 12 || in.h > 12) throw std::invalid_argument("dumpdp prints whole tables, so the image must be at most 12x12");
        e = compute_energy(to_gray(in), EnergyKind::L1);
    }
    std::printf("%s", format_dp_tables(e, dp_tables(e)).c_str());
    return 0;
}

int cmd_order(int argc, char** argv) {
    Args a(argc, argv, {"--min-width", "-o"}, {});
    Image in = load_image(a.input());
    if (!a.opt.count("--min-width")) throw UsageError("order needs --min-width K");
    int minW = arg([&] { return parse_size(a.opt["--min-width"], in.w); });
    std::string path = a.get("-o", "seams.json");
    write_seams_json(path, in.w, in.h, seams_as_cut(in.h, in.w, carve(in, minW, in.h, CarveOptions{}).removed_v));
    std::printf("seams for %dx%d down to width %d -> %s\n", in.w, in.h, minW, path.c_str());
    return 0;
}

int cmd_bench(int argc, char** argv) {
    Args a(argc, argv, {"--suite", "--out"}, {"--quick"});
    if (!a.pos.empty()) throw UsageError("bench takes no positional arguments");
    BenchArgs b;
    b.suite = a.get("--suite", "all");
    b.out = a.get("--out", "results");
    b.quick = a.flags.count("--quick") > 0;
    const std::set<std::string> suites = {"all", "scaling", "carve", "flow", "improve", "sobel", "count"};
    if (!suites.count(b.suite)) throw UsageError("unknown suite '" + b.suite + "'");
    return run_bench(b);
}

// One image of `seamcarve demo`: the files tools/make_viz.py embeds in viz/index.html. Outputs follow
// the input's format: a photo (.jpg) gives JPEGs, a drawing (.png) gives PNGs, so flat colours stay
// exact and photos stay small.
void demo_image(const std::string& file, const std::string& out, std::mutex& log) {
    Timer total;
    const std::filesystem::path path(file);
    std::string ext = path.extension().string();
    for (char& ch : ext) ch = char(std::tolower(static_cast<unsigned char>(ch)));
    const bool photo = ext == ".jpg" || ext == ".jpeg";
    const std::string dir = out + "/" + path.stem().string();
    auto save = [&](const std::string& name, const auto& img) {
        if (photo) save_jpg(dir + "/" + name + ".jpg", img);
        else save_png(dir + "/" + name + ".png", img);
    };
    Image full = load_image(file);

    // Every method on the same input (width <= 320, area-averaged), so removed energy and time are comparable.
    Image in = full.w <= 320 ? full : resize_area(full, 320, std::max(1, int(std::lround(double(full.h) * 320 / full.w))));
    save("input", in);
    Csv summary(dir + "/summary.csv", {"method", "input_w", "input_h", "out_w", "out_h", "removed_energy", "ms"});
    for (const std::string m : {"dp", "greedy", "dinic"}) {
        CarveOptions o;
        o.method = m;
        CarveResult r = carve(in, in.w / 2, in.h, o);
        save(m, r.out);
        summary.row(m, in.w, in.h, r.out.w, r.out.h, r.removed_energy, r.ms);
    }
    for (const std::string& k : baseline_names()) {
        Timer t;
        Image b = run_baseline(k, in, in.w / 2, in.h, EnergyKind::L1);
        double ms = t.ms();
        save(k, b);
        summary.row(k, in.w, in.h, b.w, b.h, std::string(""), ms);
    }

    // The slider works on the picture itself: the file, unchanged, and the seams of a DP carve to 50% width.
    // remove first: MinGW's copy_file fails with "File exists" even with overwrite_existing
    std::filesystem::remove(dir + "/original" + ext);
    std::filesystem::copy_file(file, dir + "/original" + ext);
    CarveResult r = carve(full, std::max(1, full.w / 2), full.h, CarveOptions{});
    write_seams_json(dir + "/seams.json", full.w, full.h, seams_as_cut(full.h, full.w, r.removed_v));
    save("energy", energy_view(compute_energy(to_gray(full), EnergyKind::L1)));

    std::lock_guard<std::mutex> lock(log);
    std::printf("[demo] %s done in %.1f s\n", dir.c_str(), total.ms() / 1000);
    std::fflush(stdout);
}

int cmd_demo(int argc, char** argv) {
    Args a(argc, argv, {"--out"}, {});
    if (!a.pos.empty()) throw UsageError("demo takes no positional arguments");
    std::string out = a.get("--out", "results/demo");
    std::vector<std::string> files = bench::image_files();
    if (files.empty()) throw std::runtime_error("no images in images/");
    // One thread per image: graph cut on the 320 px input takes about two minutes a photo, so the
    // images run side by side. The outputs do not depend on timing; the ms column is indicative only.
    std::mutex log;
    std::vector<std::exception_ptr> errors(files.size());
    std::vector<std::thread> workers;
    for (size_t i = 0; i < files.size(); ++i)
        workers.emplace_back([&, i] {
            try {
                demo_image(files[i], out, log);
            } catch (...) {
                errors[i] = std::current_exception();
            }
        });
    for (std::thread& t : workers) t.join();
    for (const std::exception_ptr& e : errors)
        if (e) std::rethrow_exception(e);
    return 0;
}

int cmd_methods(int, char**) {
    for (const std::string& m : finder_names()) std::printf("%s\n", m.c_str());
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    const std::map<std::string, int (*)(int, char**)> commands = {
        {"carve", cmd_carve}, {"baseline", cmd_baseline}, {"energy", cmd_energy}, {"dumpdp", cmd_dumpdp},
        {"order", cmd_order}, {"bench", cmd_bench},       {"demo", cmd_demo},     {"methods", cmd_methods}};
    PauseOnExit pause_on_exit;  // a double-clicked window must not close before the text can be read
    if (argc < 2) {
        if (opened_by_double_click())
            std::fprintf(stderr, "seamcarve is a command-line tool, so double-clicking it does nothing.\n"
                                 "To rebuild everything, double-click build_all.bat (or run: mingw32-make results).\n"
                                 "To see every command, open a terminal here (type cmd in the Explorer address bar) and run: seamcarve --help\n");
        else
            std::fprintf(stderr, "%s", kUsage);
        return 2;
    }
    std::string cmd = argv[1];
    if (cmd == "--help" || cmd == "-h" || cmd == "help") {
        std::printf("%s", kUsage);
        return 0;
    }
    auto it = commands.find(cmd);
    try {
        if (it == commands.end()) throw UsageError("unknown command '" + cmd + "'");
        return it->second(argc, argv);
    } catch (const UsageError& e) {
        std::fprintf(stderr, "error: %s\nrun 'seamcarve --help' for usage\n", e.what());
        return 2;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
