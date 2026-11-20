// energy.h - the energy map (brief 1a). Owner: Shahana.
//
// Definitions (borders are clamped: I(r, c) = g.at(clamp(r, 0, h-1), clamp(c, 0, w-1))):
//   Gx = [I(r-1,c+1) + 2I(r,c+1) + I(r+1,c+1)] - [I(r-1,c-1) + 2I(r,c-1) + I(r+1,c-1)]
//   Gy = [I(r+1,c-1) + 2I(r+1,c) + I(r+1,c+1)] - [I(r-1,c-1) + 2I(r-1,c) + I(r-1,c+1)]
//   L1 = |Gx| + |Gy| (at most 2040)    L2 = round(sqrt(Gx^2 + Gy^2))
//   Gradient = |I(r,c+1) - I(r,c-1)|   (the brief's horizontal central difference; this is NOT
//   "forward energy", which in the literature means the energy a removal inserts)
#pragma once

#include "core.h"

namespace sc {

enum class EnergyKind { L1, L2, Gradient };

inline EnergyKind parse_energy(const std::string& s) {
    if (s == "l1") return EnergyKind::L1;
    if (s == "l2") return EnergyKind::L2;
    if (s == "grad") return EnergyKind::Gradient;
    throw std::invalid_argument("unknown energy '" + s + "' (expected l1, l2 or grad)");
}

inline std::string to_string(EnergyKind k) {
    switch (k) {
        case EnergyKind::L1: return "l1";
        case EnergyKind::L2: return "l2";
        case EnergyKind::Gradient: return "grad";
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

inline Energy gradient_energy(const Gray& g) {
    Energy out(g.h, g.w);
    detail::direct_rows(g, EnergyKind::Gradient, out, 0, g.h);
    return out;
}

inline Energy compute_energy(const Gray& g, EnergyKind k) {
    Energy out(g.h, g.w);
    detail::direct_rows(g, k, out, 0, g.h);
    return out;
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
    SC_CHECK_GRID_EQ(compute_energy(g, EnergyKind::Gradient), energy_rows(1, {10, 20, 10}));
    // A single column: every horizontal difference clamps to the pixel itself.
    Gray col(3, 1);
    col.a = {0, 100, 200};
    SC_CHECK_GRID_EQ(compute_energy(col, EnergyKind::Gradient), Energy(3, 1, 0));
    // Gx = 0; Gy = 4 I(r+1) - 4 I(r-1) with clamping = [400, 800, 400].
    SC_CHECK_GRID_EQ(compute_energy(col, EnergyKind::L1), transpose(energy_rows(1, {400, 800, 400})));
    SC_CHECK_GRID_EQ(compute_energy(Gray(1, 1, 9), EnergyKind::L2), Energy(1, 1, 0));
}

}  // namespace sc
#endif
