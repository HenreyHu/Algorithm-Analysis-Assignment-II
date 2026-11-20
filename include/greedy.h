// greedy.h - greedy seam finder and the visual-quality baselines (crop, scale, column
// removal, per-row minimum removal) that seam carving is compared against.
//
// Greedy: start at the cheapest pixel of the top row, then in every row step to the cheapest of the
// (up to) three pixels below. Theta(h) per seam once the energy map exists, versus Theta(hw) for DP.
// Height reductions reuse the width code through transpose.
#pragma once

#include "core.h"
#include "energy.h"

namespace sc {

inline Seam greedy_seam(const Energy& e, Stats& st) {
    if (e.h <= 0 || e.w <= 0) throw std::invalid_argument("greedy_seam: empty energy map");
    Seam s(e.h);
    const int32_t* top = &e.at(0, 0);
    s[0] = int(std::min_element(top, top + e.w) - top);  // leftmost minimum
    for (int r = 1; r < e.h; ++r) {
        int c = s[r - 1], best = std::max(c - 1, 0);
        for (int j = best + 1; j <= std::min(c + 1, e.w - 1); ++j)
            if (e.at(r, j) < e.at(r, best)) best = j;  // order -1, 0, +1; strictly smaller wins
        s[r] = best;
    }
    st.cells += 3LL * e.h;
    return s;
}

namespace detail {

inline void check_target(const Image& in, int W, int H, const char* who) {
    if (W < 1 || H < 1 || W > in.w || H > in.h)
        throw std::invalid_argument(std::string(who) + ": target " + std::to_string(W) + "x" + std::to_string(H) +
                                    " must be between 1x1 and the input " + std::to_string(in.w) + "x" +
                                    std::to_string(in.h) + " (enlarging is out of scope)");
}

// Applies a width-only reduction to the width, then to the height through transpose.
template <class F>
Image reduce_both(const Image& in, int W, int H, F reduce_width) {
    Image out = reduce_width(in, W);
    if (H < out.h) out = transpose(reduce_width(transpose(out), H));
    return out;
}

inline Image remove_low_columns_w(const Image& in, int W, EnergyKind k) {
    if (W == in.w) return in;
    Energy e = compute_energy(to_gray(in), k);
    std::vector<int64_t> sum(in.w, 0);
    for (int r = 0; r < e.h; ++r)
        for (int c = 0; c < e.w; ++c) sum[c] += e.at(r, c);
    std::vector<int> order(in.w);
    for (int c = 0; c < in.w; ++c) order[c] = c;
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return sum[a] < sum[b]; });  // ties: leftmost
    std::vector<bool> drop(in.w, false);
    for (int i = 0; i < in.w - W; ++i) drop[order[i]] = true;
    Image out(in.h, W);
    for (int r = 0; r < in.h; ++r)
        for (int c = 0, o = 0; c < in.w; ++c)
            if (!drop[c]) out.at(r, o++) = in.at(r, c);
    return out;
}

inline Image remove_row_minima_w(Image img, int W, EnergyKind k) {
    while (img.w > W) {
        Energy e = compute_energy(to_gray(img), k);
        Image out(img.h, img.w - 1);
        for (int r = 0; r < img.h; ++r) {
            const int32_t* row = &e.at(r, 0);
            int m = int(std::min_element(row, row + img.w) - row);
            for (int c = 0, o = 0; c < img.w; ++c)
                if (c != m) out.at(r, o++) = img.at(r, c);
        }
        img = std::move(out);
    }
    return img;
}

}  // namespace detail

inline Image crop_center(const Image& in, int W, int H) {
    detail::check_target(in, W, H, "crop_center");
    const int x0 = (in.w - W) / 2, y0 = (in.h - H) / 2;
    Image out(H, W);
    for (int r = 0; r < H; ++r)
        for (int c = 0; c < W; ++c) out.at(r, c) = in.at(y0 + r, x0 + c);
    return out;
}

// Pixel centres are aligned: source x = (X + 0.5) * w / W - 0.5 (clamped), the same for y.
// When W == w and H == h every sample lands exactly on a pixel, so the input is returned unchanged.
inline Image scale_bilinear(const Image& in, int W, int H) {
    detail::check_target(in, W, H, "scale_bilinear");
    auto source = [](int X, int n_out, int n_in, int& i0, int& i1, double& f) {
        double x = std::clamp((X + 0.5) * n_in / n_out - 0.5, 0.0, double(n_in - 1));
        i0 = int(x);
        i1 = std::min(i0 + 1, n_in - 1);
        f = x - i0;
    };
    Image out(H, W);
    for (int Y = 0; Y < H; ++Y) {
        int y0, y1;
        double fy;
        source(Y, H, in.h, y0, y1, fy);
        for (int X = 0; X < W; ++X) {
            int x0, x1;
            double fx;
            source(X, W, in.w, x0, x1, fx);
            auto mix = [&](uint8_t RGB::*ch) {
                double top = (1 - fx) * (in.at(y0, x0).*ch) + fx * (in.at(y0, x1).*ch);
                double bottom = (1 - fx) * (in.at(y1, x0).*ch) + fx * (in.at(y1, x1).*ch);
                return uint8_t(std::lround((1 - fy) * top + fy * bottom));
            };
            out.at(Y, X) = {mix(&RGB::r), mix(&RGB::g), mix(&RGB::b)};
        }
    }
    return out;
}

inline Image remove_low_columns(const Image& in, int W, int H, EnergyKind k) {
    detail::check_target(in, W, H, "remove_low_columns");
    return detail::reduce_both(in, W, H, [k](const Image& img, int n) { return detail::remove_low_columns_w(img, n, k); });
}


// minimum pixel, independently of the other rows.
inline Image remove_row_minima(const Image& in, int W, int H, EnergyKind k) {
    detail::check_target(in, W, H, "remove_row_minima");
    return detail::reduce_both(in, W, H, [k](const Image& img, int n) { return detail::remove_row_minima_w(img, n, k); });
}

inline Image run_baseline(const std::string& kind, const Image& in, int W, int H, EnergyKind k) {
    if (kind == "crop") return crop_center(in, W, H);
    if (kind == "scale") return scale_bilinear(in, W, H);
    if (kind == "column") return remove_low_columns(in, W, H, k);
    if (kind == "rowmin") return remove_row_minima(in, W, H, k);
    throw std::invalid_argument("unknown baseline '" + kind + "' (expected crop, scale, column or rowmin)");
}

inline std::vector<std::string> baseline_names() { return {"crop", "scale", "column", "rowmin"}; }

}  

// ---------------------------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------------------------
#ifdef SC_TESTS
#include "dp.h"

namespace sc {

SC_TEST(GR1, "greedy seams are valid (1,000 random grids up to 30x30)") {
    std::mt19937 rng(401);
    std::uniform_int_distribution<int> dim(1, 30);
    for (int t = 0; t < 1000; ++t) {
        Energy e = random_energy(rng, dim(rng), dim(rng), t % 2 ? 3 : 2040);
        Stats st;
        std::string why;
        SC_CHECK(seam_valid(greedy_seam(e, st), e.h, e.w, &why));
    }
}

SC_TEST(GR2, "greedy cost >= DP cost (2,000 random grids): never better than optimal") {
    std::mt19937 rng(402);
    std::uniform_int_distribution<int> dim(1, 20);
    int strictly_worse = 0;
    for (int t = 0; t < 2000; ++t) {
        Energy e = random_energy(rng, dim(rng), dim(rng), t % 2 ? 3 : 2040);
        Stats st;
        int64_t g = seam_cost(e, greedy_seam(e, st)), d = seam_cost(e, dp_seam(e, st));
        SC_CHECK(g >= d);
        strictly_worse += g > d;
    }
    SC_CHECK(strictly_worse > 0);  // the comparison is not vacuous
}

SC_TEST(GR3, "counterexample [[9,1,0],[0,9,5],[0,9,9]]: greedy [2,2,1] costs 14, DP [1,0,0] costs 1") {
    Energy e = energy_from({{9, 1, 0}, {0, 9, 5}, {0, 9, 9}});
    Stats st;
    Seam g = greedy_seam(e, st), d = dp_seam(e, st);
    SC_CHECK_EQ(g, (Seam{2, 2, 1}));
    SC_CHECK_EQ(seam_cost(e, g), int64_t(14));
    SC_CHECK_EQ(d, (Seam{1, 0, 0}));
    SC_CHECK_EQ(seam_cost(e, d), int64_t(1));
}

SC_TEST(GR4, "h=1: greedy == DP (200 random); w=1: all zeros") {
    std::mt19937 rng(404);
    std::uniform_int_distribution<int> wd(1, 40);
    for (int t = 0; t < 200; ++t) {
        Energy e = random_energy(rng, 1, wd(rng), t % 2 ? 3 : 2040);
        Stats st;
        SC_CHECK_EQ(greedy_seam(e, st), dp_seam(e, st));
    }
    Stats st;
    SC_CHECK_EQ(greedy_seam(random_energy(rng, 9, 1, 2040), st), Seam(9, 0));
}

SC_TEST(GR5, "baselines return exactly W x H; crop offset; scale identity") {
    std::mt19937 rng(405);
    std::uniform_int_distribution<int> dim(1, 24);
    for (int t = 0; t < 50; ++t) {
        int h = dim(rng), w = dim(rng);
        int H = std::uniform_int_distribution<int>(1, h)(rng), W = std::uniform_int_distribution<int>(1, w)(rng);
        Image in = random_image(rng, h, w);
        for (const std::string& kind : baseline_names()) {
            Image out = run_baseline(kind, in, W, H, EnergyKind::L1);
            SC_CHECK_EQ(out.w, W);
            SC_CHECK_EQ(out.h, H);
        }
        SC_CHECK_EQ(crop_center(in, W, H).at(0, 0), in.at((h - H) / 2, (w - W) / 2));
        SC_CHECK_GRID_EQ(scale_bilinear(in, w, h), in);
    }
    bool threw = false;
    try {
        crop_center(Image(4, 4), 5, 4);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    SC_CHECK(threw);
}

SC_TEST(GR6, "column removal drops the flat column; row-min removal deletes one pixel per row") {
    // Columns are vertically constant with values [0, 100, 200, 100], so column 2 has
    // Gx = I(3) - I(1) = 0 and zero energy, while columns 0, 1, 3 have 400, 800, 400 per pixel.
    Image in(3, 4);
    const uint8_t v[4] = {0, 100, 200, 100};
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 4; ++c) in.at(r, c) = {v[c], v[c], v[c]};
    Image out = remove_low_columns(in, 3, 3, EnergyKind::L1);
    SC_CHECK_EQ(out.w, 3);
    for (int r = 0; r < 3; ++r) {
        SC_CHECK_EQ(out.at(r, 0), in.at(r, 0));
        SC_CHECK_EQ(out.at(r, 1), in.at(r, 1));
        SC_CHECK_EQ(out.at(r, 2), in.at(r, 3));
    }
    std::mt19937 rng(406);
    Image img = random_image(rng, 7, 9);
    Image rm = remove_row_minima(img, 8, 7, EnergyKind::L1);
    SC_CHECK_EQ(rm.w, 8);
    for (int r = 0; r < 7; ++r) {
        int skipped = 0;
        for (int c = 0, o = 0; c < 9; ++c) {
            if (o < 8 && rm.at(r, o) == img.at(r, c))
                ++o;
            else
                ++skipped;
        }
        SC_CHECK_EQ(skipped, 1);
    }
}

}  
#endif
