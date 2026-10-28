// core.h - shared types, image I/O, seam utilities, timing, CSV and the test framework.
// Owner: P1. Every other module includes this file and nothing else from the project.
//
// Conventions used across the project:
//   - A Grid is indexed (row, column) = (r, c) and has height h and width w.
//   - A vertical seam is a Seam s with s.size() == h, 0 <= s[r] < w and |s[r] - s[r-1]| <= 1.
//   - Energies are non-negative int32; sums along a seam use int64.
#pragma once

#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

// stb: the implementation is compiled here. Each executable has exactly one translation unit,
// so defining the implementation macros in a header is safe.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-function"
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"
#pragma GCC diagnostic pop

namespace sc {

// ---------------------------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------------------------

struct RGB {
    uint8_t r = 0, g = 0, b = 0;
    bool operator==(const RGB& o) const { return r == o.r && g == o.g && b == o.b; }
    bool operator!=(const RGB& o) const { return !(*this == o); }
};

template <class T>
struct Grid {
    int h = 0, w = 0;
    std::vector<T> a;  // row-major, a[r * w + c]

    Grid() = default;
    Grid(int h_, int w_, T v = T{}) : h(h_), w(w_), a(size_t(h_) * size_t(w_), v) {}

    T& at(int r, int c) { return a[size_t(r) * w + c]; }
    const T& at(int r, int c) const { return a[size_t(r) * w + c]; }
    bool operator==(const Grid& o) const { return h == o.h && w == o.w && a == o.a; }
    bool operator!=(const Grid& o) const { return !(*this == o); }
};

using Image = Grid<RGB>;
using Gray = Grid<uint8_t>;
using Energy = Grid<int32_t>;
using IndexGrid = Grid<int32_t>;
using Seam = std::vector<int>;

// Work counters filled by the seam finders (whichever apply) and reported by the benchmarks.
struct Stats {
    long long cells = 0;       // DP cells / graph nodes touched
    long long augment = 0;     // augmenting paths (max-flow)
    long long bfs_visits = 0;  // nodes dequeued by BFS / DFS (max-flow)
    long long phases = 0;      // Dinic phases
    double ms = 0;             // wall time

    void add(const Stats& o) {
        cells += o.cells;
        augment += o.augment;
        bfs_visits += o.bfs_visits;
        phases += o.phases;
        ms += o.ms;
    }
};

using Finder = std::function<Seam(const Energy&, Stats&)>;

constexpr int64_t INF = int64_t(1) << 50;

// ---------------------------------------------------------------------------------------------
// Printing (used by SC_CHECK_EQ; all overloads live in sc so lookup from tests always finds them)
// ---------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& os, const RGB& p) {
    return os << '(' << +p.r << ',' << +p.g << ',' << +p.b << ')';
}

// Declared before print_value: argument-dependent lookup for std::vector only searches std,
// so print_value finds these overloads only if they are already visible at its definition.
template <class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v);
template <class T>
std::ostream& operator<<(std::ostream& os, const Grid<T>& g);

template <class T>
void print_value(std::ostream& os, const T& v) {
    if constexpr (std::is_same_v<T, int8_t> || std::is_same_v<T, uint8_t>)
        os << +v;  // print as a number, not as a character
    else
        os << v;
}

template <class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
    os << '[';
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) os << ',';
        print_value(os, v[i]);
    }
    return os << ']';
}

template <class T>
std::ostream& operator<<(std::ostream& os, const Grid<T>& g) {
    return os << "Grid h=" << g.h << " w=" << g.w;
}

// Empty if equal, otherwise the dimensions or the first differing cell.
template <class T>
std::string grid_diff(const Grid<T>& x, const Grid<T>& y) {
    std::ostringstream o;
    if (x.h != y.h || x.w != y.w) {
        o << "dims differ: h=" << x.h << " w=" << x.w << " vs h=" << y.h << " w=" << y.w;
        return o.str();
    }
    for (int r = 0; r < x.h; ++r)
        for (int c = 0; c < x.w; ++c)
            if (!(x.at(r, c) == y.at(r, c))) {
                o << "h=" << x.h << " w=" << x.w << ", first difference at (" << r << ", " << c << "): ";
                print_value(o, x.at(r, c));
                o << " vs ";
                print_value(o, y.at(r, c));
                return o.str();
            }
    return {};
}

// ---------------------------------------------------------------------------------------------
// Image I/O
// ---------------------------------------------------------------------------------------------

inline void ensure_parent_dir(const std::string& path) {
    std::filesystem::path parent = std::filesystem::path(path).parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
}

inline Image load_image(const std::string& path) {
    int w = 0, h = 0, n = 0;
    unsigned char* data = stbi_load(path.c_str(), &w, &h, &n, 3);
    if (!data) throw std::runtime_error("cannot read image: " + path);
    Image img(h, w);
    for (size_t i = 0; i < img.a.size(); ++i) img.a[i] = {data[3 * i], data[3 * i + 1], data[3 * i + 2]};
    stbi_image_free(data);
    return img;
}

inline std::vector<unsigned char> rgb_bytes(const Image& img) {
    std::vector<unsigned char> buf(img.a.size() * 3);
    for (size_t i = 0; i < img.a.size(); ++i) {
        buf[3 * i] = img.a[i].r;
        buf[3 * i + 1] = img.a[i].g;
        buf[3 * i + 2] = img.a[i].b;
    }
    return buf;
}

inline void save_png(const std::string& path, const Image& img) {
    ensure_parent_dir(path);
    if (!stbi_write_png(path.c_str(), img.w, img.h, 3, rgb_bytes(img).data(), img.w * 3))
        throw std::runtime_error("cannot write image: " + path);
}

inline void save_png(const std::string& path, const Gray& img) {
    ensure_parent_dir(path);
    if (!stbi_write_png(path.c_str(), img.w, img.h, 1, img.a.data(), img.w))
        throw std::runtime_error("cannot write image: " + path);
}

// JPEG for photos (the demo's outputs of a .jpg input): a fraction of the PNG size, with no visible loss at 92.
inline void save_jpg(const std::string& path, const Image& img, int quality = 92) {
    ensure_parent_dir(path);
    if (!stbi_write_jpg(path.c_str(), img.w, img.h, 3, rgb_bytes(img).data(), quality))
        throw std::runtime_error("cannot write image: " + path);
}

inline void save_jpg(const std::string& path, const Gray& img, int quality = 92) {
    ensure_parent_dir(path);
    if (!stbi_write_jpg(path.c_str(), img.w, img.h, 1, img.a.data(), quality))
        throw std::runtime_error("cannot write image: " + path);
}

// ITU-R BT.601 luma in integer arithmetic, rounded, so every platform gives the same bytes.
inline Gray to_gray(const Image& img) {
    Gray g(img.h, img.w);
    for (size_t i = 0; i < img.a.size(); ++i) {
        const RGB& p = img.a[i];
        g.a[i] = uint8_t((299 * p.r + 587 * p.g + 114 * p.b + 500) / 1000);
    }
    return g;
}

inline Image resize_nearest(const Image& img, int W, int H) {
    if (W <= 0 || H <= 0) throw std::invalid_argument("resize_nearest: target size must be positive");
    Image out(H, W);
    for (int r = 0; r < H; ++r) {
        int sr = int(int64_t(r) * img.h / H);
        for (int c = 0; c < W; ++c) out.at(r, c) = img.at(sr, int(int64_t(c) * img.w / W));
    }
    return out;
}

// Area-average downscale: each output pixel is the mean of the source area it covers, so a photo
// shrinks without the jagged edges of nearest-neighbour. On a common integer axis source pixel i
// covers [i*to, (i+1)*to) and output pixel o covers [o*from, (o+1)*from), so every weight is an exact
// integer overlap and every platform gives the same bytes.
inline Image resize_area(const Image& img, int W, int H) {
    if (W <= 0 || H <= 0 || W > img.w || H > img.h)
        throw std::invalid_argument("resize_area: target must be between 1x1 and the input size");
    auto spans = [](int from, int to) {
        std::vector<std::vector<std::pair<int, int64_t>>> s(to);
        for (int o = 0; o < to; ++o) {
            const int64_t a = int64_t(o) * from, b = a + from;
            for (int64_t i = a / to; i * to < b; ++i) s[o].push_back({int(i), std::min(b, (i + 1) * to) - std::max(a, i * to)});
        }
        return s;
    };
    const auto sx = spans(img.w, W), sy = spans(img.h, H);
    std::vector<int64_t> rows(size_t(img.h) * W * 3);  // horizontal sums, weight img.w per pixel
    for (int r = 0; r < img.h; ++r)
        for (int c = 0; c < W; ++c)
            for (const auto& [i, wt] : sx[c]) {
                const RGB& p = img.at(r, i);
                int64_t* q = &rows[(size_t(r) * W + c) * 3];
                q[0] += p.r * wt, q[1] += p.g * wt, q[2] += p.b * wt;
            }
    const int64_t total = int64_t(img.w) * img.h;
    Image out(H, W);
    for (int r = 0; r < H; ++r)
        for (int c = 0; c < W; ++c) {
            int64_t s[3] = {0, 0, 0};
            for (const auto& [i, wt] : sy[r])
                for (int ch = 0; ch < 3; ++ch) s[ch] += rows[(size_t(i) * W + c) * 3 + ch] * wt;
            out.at(r, c) = {uint8_t((s[0] + total / 2) / total), uint8_t((s[1] + total / 2) / total), uint8_t((s[2] + total / 2) / total)};
        }
    return out;
}

// ---------------------------------------------------------------------------------------------
// Grid and seam utilities
// ---------------------------------------------------------------------------------------------

template <class T>
Grid<T> transpose(const Grid<T>& g) {
    Grid<T> out(g.w, g.h);
    for (int r = 0; r < g.h; ++r)
        for (int c = 0; c < g.w; ++c) out.at(c, r) = g.at(r, c);
    return out;
}

inline bool seam_valid(const Seam& s, int h, int w, std::string* why = nullptr) {
    auto fail = [&](const std::string& msg) {
        if (why) *why = msg;
        return false;
    };
    if (int(s.size()) != h)
        return fail("seam length " + std::to_string(s.size()) + " != height " + std::to_string(h));
    for (int r = 0; r < h; ++r) {
        if (s[r] < 0 || s[r] >= w)
            return fail("row " + std::to_string(r) + ": column " + std::to_string(s[r]) + " outside [0, " +
                        std::to_string(w) + ")");
        if (r > 0 && std::abs(s[r] - s[r - 1]) > 1)
            return fail("row " + std::to_string(r) + ": jump from " + std::to_string(s[r - 1]) + " to " +
                        std::to_string(s[r]));
    }
    return true;
}

// Copy of g with one pixel per row removed (the seam), so the width drops by one.
template <class T>
Grid<T> remove_seam(const Grid<T>& g, const Seam& s) {
    std::string why;
    if (!seam_valid(s, g.h, g.w, &why)) throw std::invalid_argument("remove_seam: invalid seam: " + why);
    Grid<T> out(g.h, g.w - 1);
    for (int r = 0; r < g.h; ++r) {
        const T* src = &g.a[size_t(r) * g.w];
        T* dst = out.a.data() + size_t(r) * out.w;
        std::copy(src, src + s[r], dst);
        std::copy(src + s[r] + 1, src + g.w, dst + s[r]);
    }
    return out;
}

inline int64_t seam_cost(const Energy& e, const Seam& s) {
    int64_t sum = 0;
    for (int r = 0; r < int(s.size()); ++r) sum += e.at(r, s[r]);
    return sum;
}

inline Energy random_energy(std::mt19937& g, int h, int w, int maxv) {
    std::uniform_int_distribution<int> d(0, maxv);
    Energy e(h, w);
    for (auto& v : e.a) v = d(g);
    return e;
}

// Uniform random start column, then a random step in {-1, 0, +1} per row (staying in range).
inline Seam random_seam(std::mt19937& g, int h, int w) {
    Seam s(h);
    if (h <= 0) return s;
    s[0] = std::uniform_int_distribution<int>(0, w - 1)(g);
    std::uniform_int_distribution<int> step(-1, 1);
    for (int r = 1; r < h; ++r) s[r] = std::clamp(s[r - 1] + step(g), 0, w - 1);
    return s;
}

inline Image random_image(std::mt19937& g, int h, int w) {
    std::uniform_int_distribution<int> d(0, 255);
    Image img(h, w);
    for (auto& p : img.a) p = {uint8_t(d(g)), uint8_t(d(g)), uint8_t(d(g))};
    return img;
}

// ---------------------------------------------------------------------------------------------
// Timing and CSV
// ---------------------------------------------------------------------------------------------

struct Timer {
    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    void start() { t0 = std::chrono::steady_clock::now(); }
    double ms() const { return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count(); }
};

inline double median(std::vector<double> v) {
    if (v.empty()) return 0;
    std::sort(v.begin(), v.end());
    size_t n = v.size();
    return n % 2 ? v[n / 2] : (v[n / 2 - 1] + v[n / 2]) / 2;
}

// Comma-separated file written with <cstdio>; a field with a comma, quote or newline is quoted.
class Csv {
public:
    // `note`, when given, becomes a first "# ..." line (for example the machine the numbers came from).
    Csv(const std::string& path, const std::vector<std::string>& header, const std::string& note = "") : path_(path) {
        ensure_parent_dir(path);
        f_ = std::fopen(path.c_str(), "w");
        if (!f_) throw std::runtime_error("cannot write csv: " + path);
        if (!note.empty()) std::fprintf(f_, "# %s\n", note.c_str());
        for (size_t i = 0; i < header.size(); ++i) std::fprintf(f_, "%s%s", i ? "," : "", field(header[i]).c_str());
        std::fputc('\n', f_);
    }
    Csv(const Csv&) = delete;
    Csv& operator=(const Csv&) = delete;
    ~Csv() {
        if (f_) std::fclose(f_);
    }

    template <class... A>
    void row(const A&... cols) {
        bool first = true;
        ((std::fprintf(f_, "%s%s", first ? "" : ",", field(cols).c_str()), first = false), ...);
        std::fputc('\n', f_);
        std::fflush(f_);
    }

private:
    static std::string field(const std::string& s) {
        if (s.find_first_of(",\"\n") == std::string::npos) return s;
        std::string q = "\"";
        for (char ch : s) q += ch == '"' ? std::string("\"\"") : std::string(1, ch);
        return q + "\"";
    }
    static std::string field(const char* s) { return field(std::string(s)); }
    static std::string field(double v) {
        char buf[64];
        std::snprintf(buf, sizeof buf, "%.6g", v);
        return buf;
    }
    template <class T, class = std::enable_if_t<std::is_integral_v<T>>>
    static std::string field(T v) {
        return std::to_string(v);
    }

    std::string path_;
    std::FILE* f_ = nullptr;
};

// ---------------------------------------------------------------------------------------------
// Test framework: SC_TEST self-registers, so no central list has to be edited.
// ---------------------------------------------------------------------------------------------

struct TestCase {
    std::string id, name;
    void (*fn)();
};

inline std::vector<TestCase>& test_registry() {
    static std::vector<TestCase> r;
    return r;
}

inline bool register_test(TestCase t) {
    test_registry().push_back(std::move(t));
    return true;
}

struct CheckFailure : std::runtime_error {
    using std::runtime_error::runtime_error;
};

inline std::string check_where(const char* file, int line) {
    return std::filesystem::path(file).filename().string() + ":" + std::to_string(line);
}

}  // namespace sc

#define SC_TEST(ID, NAME)                                                     \
    static void ID##_fn();                                                    \
    static const bool ID##_reg = ::sc::register_test({#ID, NAME, &ID##_fn}); \
    static void ID##_fn()

#define SC_CHECK(c)                                                                                 \
    do {                                                                                            \
        if (!(c)) throw ::sc::CheckFailure(::sc::check_where(__FILE__, __LINE__) + "  CHECK(" #c ")"); \
    } while (0)

#define SC_CHECK_EQ(a, b)                                                    \
    do {                                                                     \
        const auto _sc_a = (a); /* by value: (a) may be a reference into a temporary */ \
        const auto _sc_b = (b);                                              \
        if (!(_sc_a == _sc_b)) {                                             \
            std::ostringstream _sc_o;                                        \
            _sc_o << ::sc::check_where(__FILE__, __LINE__) << "  " #a " == " #b "  got "; \
            ::sc::print_value(_sc_o, _sc_a);                                 \
            _sc_o << " vs ";                                                 \
            ::sc::print_value(_sc_o, _sc_b);                                 \
            throw ::sc::CheckFailure(_sc_o.str());                           \
        }                                                                    \
    } while (0)

#define SC_CHECK_GRID_EQ(a, b)                                                                   \
    do {                                                                                         \
        std::string _sc_d = ::sc::grid_diff((a), (b));                                           \
        if (!_sc_d.empty())                                                                      \
            throw ::sc::CheckFailure(::sc::check_where(__FILE__, __LINE__) + "  " #a " == " #b "  " + _sc_d); \
    } while (0)

// ---------------------------------------------------------------------------------------------
// Double-click support (Windows)
// ---------------------------------------------------------------------------------------------
#ifdef _WIN32
// Three kernel32 calls, declared here so that windows.h (whose RGB macro clashes with our struct) stays out.
extern "C" __declspec(dllimport) unsigned long __stdcall GetConsoleProcessList(unsigned long* list, unsigned long count);
extern "C" __declspec(dllimport) void* __stdcall GetStdHandle(unsigned long which);
extern "C" __declspec(dllimport) int __stdcall GetConsoleMode(void* handle, unsigned long* mode);
#endif

namespace sc {

// True when Explorer started the program: it then opens a console window that belongs to this program
// alone, and the window closes the moment the program exits. In a terminal the console is shared with the
// shell, and under a pipe, a file, NUL or CI there is no console to keep open, so this is false.
inline bool opened_by_double_click() {
#ifdef _WIN32
    const unsigned long kStdInput = static_cast<unsigned long>(-10);  // STD_INPUT_HANDLE
    unsigned long mode, pids[2];
    if (std::getenv("CI") || !GetConsoleMode(GetStdHandle(kStdInput), &mode)) return false;  // stdin is not a console
    return GetConsoleProcessList(pids, 2) == 1;
#else
    return false;
#endif
}

// Declare one at the top of main(): a double-clicked window then stays open until Enter is pressed.
struct PauseOnExit {
    ~PauseOnExit() {
        if (!opened_by_double_click()) return;
        std::printf("\nPress Enter to close this window...");
        std::fflush(stdout);
        std::getchar();
    }
};

}  // namespace sc

// ---------------------------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------------------------
#ifdef SC_TESTS
namespace sc {

SC_TEST(CORE1, "transpose twice is the identity (50 grids incl. 1xn, nx1)") {
    std::mt19937 rng(101);
    std::uniform_int_distribution<int> dim(1, 12);
    for (int t = 0; t < 50; ++t) {
        int h = t == 0 ? 1 : dim(rng), w = t == 1 ? 1 : dim(rng);
        if (t == 2) h = 1, w = 9;
        if (t == 3) h = 9, w = 1;
        Energy e = random_energy(rng, h, w, 1000);
        Energy tt = transpose(e);
        SC_CHECK_EQ(tt.h, w);
        SC_CHECK_EQ(tt.w, h);
        SC_CHECK_GRID_EQ(transpose(tt), e);
    }
}

SC_TEST(CORE2, "remove_seam on a 3x4 grid removes exactly the seam pixels") {
    Energy g(3, 4);
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 4; ++c) g.at(r, c) = r * 10 + c;
    Energy got = remove_seam(g, {1, 2, 3});
    Energy want(3, 3);
    want.a = {0, 2, 3, 10, 11, 13, 20, 21, 22};
    SC_CHECK_GRID_EQ(got, want);
}

SC_TEST(CORE3, "seam_valid accepts a valid seam and rejects the three kinds of bad seam") {
    std::string why;
    SC_CHECK(seam_valid({1, 2, 1}, 3, 3, &why));
    SC_CHECK(why.empty());
    for (const Seam& bad : {Seam{0, 2, 1}, Seam{0, 1, 3}, Seam{0, 1}}) {
        why.clear();
        SC_CHECK(!seam_valid(bad, 3, 3, &why));
        SC_CHECK(!why.empty());
    }
}

SC_TEST(CORE4, "seam_cost sums the energy under the seam") {
    Energy g(3, 4);
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 4; ++c) g.at(r, c) = r * 10 + c;
    SC_CHECK_EQ(seam_cost(g, {1, 2, 3}), int64_t(36));
}

SC_TEST(CORE5, "PNG write then read gives the same pixels") {
    std::mt19937 rng(105);
    Image img = random_image(rng, 5, 7);  // W x H = 7 x 5
    save_png("tests_tmp/rt.png", img);
    SC_CHECK_GRID_EQ(load_image("tests_tmp/rt.png"), img);
}

SC_TEST(CORE6, "to_gray uses rounded integer BT.601 weights") {
    Image img(1, 5);
    img.a = {{255, 0, 0}, {0, 255, 0}, {0, 0, 255}, {255, 255, 255}, {0, 0, 0}};
    Gray want(1, 5);
    want.a = {76, 150, 29, 255, 0};
    SC_CHECK_GRID_EQ(to_gray(img), want);
}

SC_TEST(CORE7, "resize_area averages exactly the covered area (hand-checked 4->2 and 3->2, identity, flat)") {
    auto row = [](std::vector<uint8_t> v) {
        Image img(1, int(v.size()));
        for (size_t i = 0; i < v.size(); ++i) img.a[i] = {v[i], 0, 0};
        return img;
    };
    SC_CHECK_GRID_EQ(resize_area(row({0, 10, 20, 40}), 2, 1), row({5, 30}));
    SC_CHECK_GRID_EQ(resize_area(row({0, 30, 90}), 2, 1), row({10, 70}));  // (0*2 + 30) / 3, (30 + 90*2) / 3
    SC_CHECK_GRID_EQ(resize_area(row({0, 1}), 1, 1), row({1}));            // a half rounds up
    std::mt19937 rng(507);
    Image any = random_image(rng, 5, 7);
    SC_CHECK_GRID_EQ(resize_area(any, 7, 5), any);
    Image flat(5, 7, RGB{12, 200, 99});
    SC_CHECK_GRID_EQ(resize_area(flat, 3, 2), Image(2, 3, RGB{12, 200, 99}));
}

}  // namespace sc
#endif
