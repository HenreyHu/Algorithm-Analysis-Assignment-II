// energy.h - the energy map (brief 1a), the divide-and-conquer evidence for Q2b and the
// energy-update improvement for 5d. Owner: Shahana.
//
// Definitions (borders are clamped: I(r, c) = g.at(clamp(r, 0, h-1), clamp(c, 0, w-1))):
//   Gx = [I(r-1,c+1) + 2I(r,c+1) + I(r+1,c+1)] - [I(r-1,c-1) + 2I(r,c-1) + I(r+1,c-1)]
//   Gy = [I(r+1,c-1) + 2I(r+1,c) + I(r+1,c+1)] - [I(r-1,c-1) + 2I(r-1,c) + I(r-1,c+1)]
//   L1 = |Gx| + |Gy| (at most 2040)    L2 = round(sqrt(Gx^2 + Gy^2))
//   Gradient = |I(r,c+1) - I(r,c-1)|   (the brief's horizontal central difference; this is NOT
//   "forward energy", which in the literature means the energy a removal inserts)
//
// Complexity and Q2b (is divide and conquer useful for the convolution?):
//   - Direct Sobel: 2 x 9 multiply-adds per pixel, Theta(hw).
//   - Separable Sobel: [1,2,1] smoothing times [-1,0,1] difference, 2 x (3 + 3) per pixel,
//     still Theta(hw) with a smaller constant. Identical output, because clamping is per axis.
//   - Recursive row-halving D&C (sobel_dc): T(h) = 2T(h/2) + Theta(w) with a Theta(16w) base,
//     so T = Theta(hw) by the Master theorem (leaf-dominated). Correct and natural, same bound.
//   - FFT convolution is the other D&C route, Theta(N log N): it wins only for large kernels.
//     With a 3x3 kernel direct is Theta(9N) < Theta(N log N), so D&C gives no asymptotic gain.
//   - Row-band threads (compute_energy with threads > 1) are parallelism, not a better algorithm:
//     the same Theta(hw) work split across cores. D&C makes the bands independent, which is what
//     lets them run in parallel.
#pragma once

#include <thread>

#include "core.h"

namespace sc {

enum class EnergyKind { L1, L2, Gradient };
enum class UpdateMode { Full, Local, Stale };

inline EnergyKind parse_energy(const std::string& s) {
    if (s == "l1") return EnergyKind::L1;
    if (s == "l2") return EnergyKind::L2;
    if (s == "grad") return EnergyKind::Gradient;
    throw std::invalid_argument("unknown energy '" + s + "' (expected l1, l2 or grad)");
}

inline UpdateMode parse_update(const std::string& s) {
    if (s == "full") return UpdateMode::Full;
    if (s == "local") return UpdateMode::Local;
    if (s == "stale") return UpdateMode::Stale;
    throw std::invalid_argument("unknown update mode '" + s + "' (expected full, local or stale)");
}

inline std::string to_string(EnergyKind k) {
    switch (k) {
        case EnergyKind::L1: return "l1";
        case EnergyKind::L2: return "l2";
        case EnergyKind::Gradient: return "grad";
    }
    return "?";
}

inline std::string to_string(UpdateMode m) {
    switch (m) {
        case UpdateMode::Full: return "full";
        case UpdateMode::Local: return "local";
        case UpdateMode::Stale: return "stale";
    }
    return "?";
}

namespace detail {

inline int px(const Gray& g, int r, int c) {
    return g.at(std::clamp(r, 0, g.h - 1), std::clamp(c, 0, g.w - 1));
}

inline int32_t combine(int gx, int gy, EnergyKind k) {
    if (k == EnergyKind::L1) return std::abs(gx) + std::abs(gy);
    return int32_t(std::lround(std::sqrt(double(gx) * gx + double(gy) * gy)));
}

// One pixel by the definition: 9 taps per kernel for L1/L2, 2 taps for Gradient.
inline int32_t energy_at(const Gray& g, EnergyKind k, int r, int c) {
    if (k == EnergyKind::Gradient) return std::abs(px(g, r, c + 1) - px(g, r, c - 1));
    int gx = (px(g, r - 1, c + 1) + 2 * px(g, r, c + 1) + px(g, r + 1, c + 1)) -
             (px(g, r - 1, c - 1) + 2 * px(g, r, c - 1) + px(g, r + 1, c - 1));
    int gy = (px(g, r + 1, c - 1) + 2 * px(g, r + 1, c) + px(g, r + 1, c + 1)) -
             (px(g, r - 1, c - 1) + 2 * px(g, r - 1, c) + px(g, r - 1, c + 1));
    return combine(gx, gy, k);
}

inline void direct_rows(const Gray& g, EnergyKind k, Energy& out, int r0, int r1) {
    for (int r = r0; r < r1; ++r)
        for (int c = 0; c < g.w; ++c) out.at(r, c) = energy_at(g, k, r, c);
}

// Separable Sobel on rows [r0, r1). Pass 1 builds, for the band rows r0-1 .. r1 (clamped),
//   D(r,c) = I(r,c+1) - I(r,c-1)   and   S(r,c) = I(r,c-1) + 2I(r,c) + I(r,c+1);
// pass 2 combines three band rows: Gx = D(r-1) + 2D(r) + D(r+1), Gy = S(r+1) - S(r-1).
inline void separable_rows(const Gray& g, EnergyKind k, Energy& out, int r0, int r1) {
    if (r0 >= r1 || g.w == 0) return;
    const int w = g.w, n = r1 - r0 + 2;
    std::vector<int32_t> D(size_t(n) * w), S(size_t(n) * w);
    for (int b = 0; b < n; ++b) {
        const uint8_t* row = &g.a[size_t(std::clamp(r0 - 1 + b, 0, g.h - 1)) * w];
        int32_t* d = &D[size_t(b) * w];
        int32_t* s = &S[size_t(b) * w];
        for (int c = 0; c < w; ++c) {
            int left = row[std::max(c - 1, 0)], mid = row[c], right = row[std::min(c + 1, w - 1)];
            d[c] = right - left;
            s[c] = left + 2 * mid + right;
        }
    }
    for (int r = r0; r < r1; ++r) {
        size_t up = size_t(r - r0) * w, me = up + w, down = me + w;
        for (int c = 0; c < w; ++c) {
            int gx = D[up + c] + 2 * D[me + c] + D[down + c];
            int gy = S[down + c] - S[up + c];
            out.at(r, c) = combine(gx, gy, k);
        }
    }
}

inline void require_sobel(EnergyKind k, const char* who) {
    if (k == EnergyKind::Gradient) throw std::invalid_argument(std::string(who) + ": needs l1 or l2");
}

}  // namespace detail

inline Energy sobel_direct(const Gray& g, EnergyKind k) {
    detail::require_sobel(k, "sobel_direct");
    Energy out(g.h, g.w);
    detail::direct_rows(g, k, out, 0, g.h);
    return out;
}

inline Energy sobel_separable(const Gray& g, EnergyKind k) {
    detail::require_sobel(k, "sobel_separable");
    Energy out(g.h, g.w);
    detail::separable_rows(g, k, out, 0, g.h);
    return out;
}

// Divide and conquer over rows: solve(r0, r1) splits at the middle row until a band has at most
// 16 rows, then computes it directly. Each half reads its one-row halo straight from the shared
// input, so combining is free. T(h) = 2T(h/2) + Theta(w) = Theta(hw).
inline Energy sobel_dc(const Gray& g, EnergyKind k) {
    detail::require_sobel(k, "sobel_dc");
    Energy out(g.h, g.w);
    std::function<void(int, int)> solve = [&](int r0, int r1) {
        if (r1 - r0 <= 16) {
            detail::direct_rows(g, k, out, r0, r1);
            return;
        }
        int mid = r0 + (r1 - r0) / 2;
        solve(r0, mid);
        solve(mid, r1);
    };
    solve(0, g.h);
    return out;
}

inline Energy gradient_energy(const Gray& g) {
    Energy out(g.h, g.w);
    detail::direct_rows(g, EnergyKind::Gradient, out, 0, g.h);
    return out;
}

// Rows [r0, r1) of the energy map (separable for L1/L2), the unit of work for threading.
inline void compute_energy_rows(const Gray& g, EnergyKind k, Energy& out, int r0, int r1) {
    if (k == EnergyKind::Gradient)
        detail::direct_rows(g, k, out, r0, r1);
    else
        detail::separable_rows(g, k, out, r0, r1);
}

inline Energy compute_energy(const Gray& g, EnergyKind k, int threads = 1) {
    Energy out(g.h, g.w);
    int t = std::clamp(threads, 1, std::max(g.h, 1));
    if (t == 1) {
        compute_energy_rows(g, k, out, 0, g.h);
        return out;
    }
    // t equal row bands; band 0 runs on the calling thread.
    std::vector<std::thread> pool;
    for (int b = 1; b < t; ++b)
        pool.emplace_back(compute_energy_rows, std::cref(g), k, std::ref(out), int(int64_t(g.h) * b / t),
                          int(int64_t(g.h) * (b + 1) / t));
    compute_energy_rows(g, k, out, 0, int(int64_t(g.h) / t));
    for (auto& th : pool) th.join();
    return out;
}

// Removes seam s from the gray image and its energy map, then brings the energy up to date:
//   Full  - recompute everything, Theta(hw).
//   Local - recompute only columns [lo-1, hi] of each row r, where lo/hi are the min/max seam
//           column over rows r-1, r, r+1: the only pixels whose 3x3 neighbourhood changed.
//           Theta(h) per seam, and equal to Full (test EN5).
//   Stale - keep the shifted old values (negative control, test EN6).
inline void remove_and_update(Gray& g, Energy& e, const Seam& s, EnergyKind k, UpdateMode m,
                              Stats* st = nullptr) {
    g = remove_seam(g, s);
    e = remove_seam(e, s);
    long long cells = 0;
    if (m == UpdateMode::Full) {
        e = compute_energy(g, k);
        cells = (long long)g.h * g.w;
    } else if (m == UpdateMode::Local) {
        for (int r = 0; r < g.h; ++r) {
            int lo = s[r], hi = s[r];
            if (r > 0) lo = std::min(lo, s[r - 1]), hi = std::max(hi, s[r - 1]);
            if (r + 1 < g.h) lo = std::min(lo, s[r + 1]), hi = std::max(hi, s[r + 1]);
            int c0 = std::max(lo - 1, 0), c1 = std::min(hi, g.w - 1);
            for (int c = c0; c <= c1; ++c) e.at(r, c) = detail::energy_at(g, k, r, c);
            cells += std::max(0, c1 - c0 + 1);
        }
    }
    if (st) st->cells += cells;
}

// Linear scale to 0..255 for viewing: 255 * e / max(e), all 0 if the map is flat.
inline Gray energy_view(const Energy& e) {
    Gray out(e.h, e.w);
    int32_t mx = e.a.empty() ? 0 : *std::max_element(e.a.begin(), e.a.end());
    if (mx > 0)
        for (size_t i = 0; i < e.a.size(); ++i) out.a[i] = uint8_t(int64_t(255) * e.a[i] / mx);
    return out;
}

}  // namespace sc

// ---------------------------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------------------------
#ifdef SC_TESTS
namespace sc {

inline Gray random_gray(std::mt19937& rng, int h, int w) {
    std::uniform_int_distribution<int> d(0, 255);
    Gray g(h, w);
    for (auto& v : g.a) v = uint8_t(d(rng));
    return g;
}

inline Gray gray_rows(int h, const std::vector<int>& row) {
    Gray g(h, int(row.size()));
    for (int r = 0; r < h; ++r)
        for (int c = 0; c < g.w; ++c) g.at(r, c) = uint8_t(row[c]);
    return g;
}

inline Energy energy_rows(int h, const std::vector<int>& row) {
    Energy e(h, int(row.size()));
    for (int r = 0; r < h; ++r)
        for (int c = 0; c < e.w; ++c) e.at(r, c) = row[c];
    return e;
}

const EnergyKind kAllKinds[] = {EnergyKind::L1, EnergyKind::L2, EnergyKind::Gradient};

SC_TEST(EN1, "constant image has zero energy (L1, L2, Gradient)") {
    for (auto [h, w] : {std::pair{1, 1}, std::pair{1, 5}, std::pair{4, 1}, std::pair{6, 7}})
        for (EnergyKind k : kAllKinds) SC_CHECK_GRID_EQ(compute_energy(Gray(h, w, 77), k), Energy(h, w, 0));
}

SC_TEST(EN2, "vertical edge [0,0,255,255]: L1 = L2 = [0,1020,1020,0], Gradient = [0,255,255,0]") {
    Gray g = gray_rows(3, {0, 0, 255, 255});
    SC_CHECK_GRID_EQ(compute_energy(g, EnergyKind::L1), energy_rows(3, {0, 1020, 1020, 0}));
    SC_CHECK_GRID_EQ(compute_energy(g, EnergyKind::L2), energy_rows(3, {0, 1020, 1020, 0}));
    SC_CHECK_GRID_EQ(compute_energy(g, EnergyKind::Gradient), energy_rows(3, {0, 255, 255, 0}));
}

SC_TEST(EN3, "h=1 and w=1 use clamped borders: [0,10,20] gives L1 [40,80,40], Gradient [10,20,10]") {
    Gray g = gray_rows(1, {0, 10, 20});
    SC_CHECK_GRID_EQ(compute_energy(g, EnergyKind::L1), energy_rows(1, {40, 80, 40}));
    SC_CHECK_GRID_EQ(sobel_direct(g, EnergyKind::L1), energy_rows(1, {40, 80, 40}));
    SC_CHECK_GRID_EQ(sobel_dc(g, EnergyKind::L1), energy_rows(1, {40, 80, 40}));
    SC_CHECK_GRID_EQ(compute_energy(g, EnergyKind::Gradient), energy_rows(1, {10, 20, 10}));
    // A single column: every horizontal difference clamps to the pixel itself.
    Gray col(3, 1);
    col.a = {0, 100, 200};
    SC_CHECK_GRID_EQ(compute_energy(col, EnergyKind::Gradient), Energy(3, 1, 0));
    // Gx = 0; S(r) = 4 I(r), so Gy = S(r+1) - S(r-1) with clamping = [400, 800, 400].
    SC_CHECK_GRID_EQ(compute_energy(col, EnergyKind::L1), transpose(energy_rows(1, {400, 800, 400})));
    SC_CHECK_GRID_EQ(compute_energy(Gray(1, 1, 9), EnergyKind::L2), Energy(1, 1, 0));
}

SC_TEST(EN4, "direct == separable == recursive D&C == naive 2-D convolution (200 random, L1/L2)") {
    const int KX[3][3] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};
    const int KY[3][3] = {{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}};
    std::mt19937 rng(204);
    std::uniform_int_distribution<int> dim(1, 40);
    for (int t = 0; t < 200; ++t) {
        Gray g = random_gray(rng, dim(rng), dim(rng));
        for (EnergyKind k : {EnergyKind::L1, EnergyKind::L2}) {
            Energy naive(g.h, g.w);
            for (int r = 0; r < g.h; ++r)
                for (int c = 0; c < g.w; ++c) {
                    int gx = 0, gy = 0;
                    for (int dr = -1; dr <= 1; ++dr)
                        for (int dc = -1; dc <= 1; ++dc) {
                            int v = g.at(std::clamp(r + dr, 0, g.h - 1), std::clamp(c + dc, 0, g.w - 1));
                            gx += KX[dr + 1][dc + 1] * v;
                            gy += KY[dr + 1][dc + 1] * v;
                        }
                    naive.at(r, c) = k == EnergyKind::L1 ? std::abs(gx) + std::abs(gy)
                                                         : int32_t(std::lround(std::hypot(double(gx), double(gy))));
                }
            SC_CHECK_GRID_EQ(sobel_direct(g, k), naive);
            SC_CHECK_GRID_EQ(sobel_separable(g, k), naive);
            SC_CHECK_GRID_EQ(sobel_dc(g, k), naive);
        }
    }
}

// Removes 5 random seams three ways (Full, Local, Stale) from copies of the same image.
// Calls check(full, local, stale) after every removal.
template <class F>
void run_update_chains(uint32_t seed, int images, F check) {
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> hd(1, 12), wd(6, 14);
    for (int t = 0; t < images; ++t) {
        Gray g0 = random_gray(rng, hd(rng), wd(rng));
        for (EnergyKind k : kAllKinds) {
            Gray gf = g0, gl = g0, gs = g0;
            Energy ef = compute_energy(g0, k), el = ef, es = ef;
            for (int step = 0; step < 5; ++step) {
                Seam s = random_seam(rng, gf.h, gf.w);
                remove_and_update(gf, ef, s, k, UpdateMode::Full);
                remove_and_update(gl, el, s, k, UpdateMode::Local);
                remove_and_update(gs, es, s, k, UpdateMode::Stale);
                check(ef, el, es);
            }
        }
    }
}

SC_TEST(EN5, "Local update == Full recompute after each of 5 removals (300 random, all kinds)") {
    run_update_chains(205, 300, [](const Energy& full, const Energy& local, const Energy&) {
        SC_CHECK_GRID_EQ(local, full);
    });
}

SC_TEST(EN6, "Stale update differs from Full on some case (the negative control is not a no-op)") {
    int differ = 0;
    run_update_chains(205, 300, [&](const Energy& full, const Energy&, const Energy& stale) {
        differ += stale != full;
    });
    SC_CHECK(differ > 0);
}

SC_TEST(EN7, "4-thread energy == 1-thread energy (incl. fewer rows than threads)") {
    std::mt19937 rng(207);
    std::uniform_int_distribution<int> hd(1, 300), wd(1, 40);
    for (int t = 0; t < 40; ++t) {
        Gray g = random_gray(rng, t < 4 ? t + 1 : hd(rng), wd(rng));
        for (EnergyKind k : kAllKinds) SC_CHECK_GRID_EQ(compute_energy(g, k, 4), compute_energy(g, k, 1));
    }
}

}  // namespace sc
#endif
