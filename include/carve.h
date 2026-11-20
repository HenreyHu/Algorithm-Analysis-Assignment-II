// carve.h - the carving driver (brief 1c): find a seam, remove it, update the energy, repeat, in
// both directions; plus the seam overlay, step frames, the removal-order map and the packed seams
// behind the visualiser's slider (bonus 7a/7b). Owner: Yan Xin.
//
// Index tracking: idx starts as idx(r, c) = r*w0 + c (the pixel's original index) and goes through
// exactly the same remove_seam/transpose calls as the image, so at any time idx tells which
// original pixel sits where. That is how every removed seam is recorded in original coordinates.
//
// Horizontal seams reuse the vertical code on the transposed image. The energy map can be
// transposed instead of recomputed for L1/L2: under transpose Sobel's Gx and Gy swap exactly, so
// |Gx| + |Gy| and sqrt(Gx^2 + Gy^2) are unchanged. The Gradient energy only looks left and right,
// so it is not symmetric and is recomputed from the transposed gray image (test CV6).
#pragma once

#include <algorithm>
#include <cstdio>
#include <numeric>

#include "core.h"
#include "dp.h"
#include "energy.h"
#include "graphcut.h"
#include "greedy.h"

namespace sc {

struct CarveOptions {
    std::string method = "dp";
    EnergyKind energy = EnergyKind::L1;
    UpdateMode update = UpdateMode::Local;
    int threads = 1;          // energy threads; "mt" uses max(threads, 2) DP threads
    std::string frames_dir;   // frames off when empty or frame_every == 0
    int frame_every = 0;
};

struct CarveResult {
    Image out;
    IndexGrid kept;                                         // original index r*w0 + c of each output pixel
    std::vector<std::vector<int32_t>> removed_v, removed_h;  // original indices of each removed seam
    std::vector<int64_t> seam_costs;
    int64_t removed_energy = 0;
    double ms = 0;
    Stats stats;
};

inline std::vector<std::string> finder_names() {
    return {"dp", "rolling", "memo", "mt", "greedy", "brute", "ff", "ek", "dinic"};
}

inline Finder get_finder(const std::string& name, int threads = 1) {
    if (name == "dp") return dp_seam;
    if (name == "rolling") return dp_seam_rolling;
    if (name == "memo") return dp_seam_memo;
    if (name == "mt") {
        int t = std::max(threads, 2);
        return [t](const Energy& e, Stats& st) { return dp_seam_mt(e, st, t); };
    }
    if (name == "greedy") return greedy_seam;
    if (name == "brute") return brute_seam;
    if (name == "ff") return ff_seam;
    if (name == "ek") return ek_seam;
    if (name == "dinic") return dinic_seam;
    std::string all;
    for (const auto& n : finder_names()) all += (all.empty() ? "" : ", ") + n;
    throw std::invalid_argument("unknown method '" + name + "' (expected one of: " + all + ")");
}

inline CarveResult carve(const Image& in, int W, int H, const CarveOptions& o) {
    if (W < 1 || H < 1 || W > in.w || H > in.h)
        throw std::invalid_argument("carve: target " + std::to_string(W) + "x" + std::to_string(H) +
                                    " must be between 1x1 and the input " + std::to_string(in.w) + "x" +
                                    std::to_string(in.h) + " (enlarging is out of scope)");
    const Finder find = get_finder(o.method, o.threads);
    if ((o.method == "ff" || o.method == "ek" || o.method == "dinic") && (long long)in.w * in.h > 40000)
        std::fprintf(stderr, "warning: graph cut on %dx%d is slow; consider resizing to <= 200x200 first\n", in.w, in.h);

    Timer total;
    CarveResult res;
    Image img = in;
    IndexGrid idx(in.h, in.w);
    for (int i = 0; i < int(idx.a.size()); ++i) idx.a[i] = i;
    Gray G = to_gray(img);
    Energy E = compute_energy(G, o.energy, o.threads);
    int step = 0;

    // Removes seams until the (current) width reaches target; `transposed` only affects frames.
    auto pass = [&](int target, bool transposed, std::vector<std::vector<int32_t>>& removed) {
        while (img.w > target) {
            Stats s;
            Timer t;
            Seam seam = find(E, s);
            s.ms = t.ms();
            res.stats.add(s);
            int64_t cost = seam_cost(E, seam);
            res.seam_costs.push_back(cost);
            res.removed_energy += cost;
            std::vector<int32_t> ids(img.h);
            for (int r = 0; r < img.h; ++r) ids[r] = idx.at(r, seam[r]);
            removed.push_back(std::move(ids));
            ++step;
            if (!o.frames_dir.empty() && o.frame_every > 0 && step % o.frame_every == 0) {
                Image frame = img;
                for (int r = 0; r < img.h; ++r) frame.at(r, seam[r]) = {255, 0, 0};
                char name[32];
                std::snprintf(name, sizeof name, "/step_%04d.png", step);
                save_png(o.frames_dir + name, transposed ? transpose(frame) : frame);
            }
            img = remove_seam(img, seam);
            idx = remove_seam(idx, seam);
            remove_and_update(G, E, seam, o.energy, o.update, &res.stats);
        }
    };

    pass(W, false, res.removed_v);
    if (H < img.h) {
        img = transpose(img);
        idx = transpose(idx);
        G = transpose(G);
        E = o.energy == EnergyKind::Gradient ? compute_energy(G, o.energy, o.threads) : transpose(E);
        pass(H, true, res.removed_h);
        img = transpose(img);
        idx = transpose(idx);
    }
    res.out = std::move(img);
    res.kept = std::move(idx);
    res.ms = total.ms();
    return res;
}

// All removed seams on the original image at once: vertical seams red, horizontal seams blue.
inline Image overlay_seams(const Image& original, const CarveResult& r) {
    Image out = original;
    for (const auto& seam : r.removed_v)
        for (int32_t i : seam) out.a[i] = {255, 0, 0};
    for (const auto& seam : r.removed_h)
        for (int32_t i : seam) out.a[i] = {0, 0, 255};
    return out;
}

// order(r, c) = k if pixel (r, c) was removed by the k-th vertical seam (0-based) while carving
// from w down to minW, or -1 if it is never removed.
inline IndexGrid removal_order(const Image& in, int minW, const CarveOptions& o) {
    CarveResult r = carve(in, minW, in.h, o);
    IndexGrid order(in.h, in.w, -1);
    for (int k = 0; k < int(r.removed_v.size()); ++k)
        for (int32_t i : r.removed_v[k]) order.a[i] = k;
    return order;
}

// The image after k vertical seams, rebuilt from the order map alone (no seam search): keep the
// pixels with order == -1 or order >= k, in their original order. The slider does exactly this.
inline Image width_from_order(const Image& in, const IndexGrid& order, int k) {
    Image out(in.h, in.w - k);
    for (int r = 0; r < in.h; ++r) {
        int o = 0;
        for (int c = 0; c < in.w; ++c) {
            int32_t v = order.at(r, c);
            if (v == -1 || v >= k) {
                if (o == out.w) throw std::logic_error("width_from_order: row " + std::to_string(r) + " keeps too many pixels");
                out.at(r, o++) = in.at(r, c);
            }
        }
        if (o != out.w) throw std::logic_error("width_from_order: row " + std::to_string(r) + " keeps too few pixels");
    }
    return out;
}

inline std::string base64(const std::vector<unsigned char>& bytes) {
    static const char* tbl = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    std::string out;
    out.reserve((bytes.size() + 2) / 3 * 4);
    for (size_t i = 0; i < bytes.size(); i += 3) {
        uint32_t v = uint32_t(bytes[i]) << 16;
        if (i + 1 < bytes.size()) v |= uint32_t(bytes[i + 1]) << 8;
        if (i + 2 < bytes.size()) v |= bytes[i + 2];
        out += tbl[(v >> 18) & 63];
        out += tbl[(v >> 12) & 63];
        out += i + 1 < bytes.size() ? tbl[(v >> 6) & 63] : '=';
        out += i + 2 < bytes.size() ? tbl[v & 63] : '=';
    }
    return out;
}

// The seams of a vertical carve of an h x w image, each in the columns of the image it was cut from
// (seam k from width w-k), recovered from the original indices in CarveResult::removed_v. A row's
// surviving columns stay sorted, so each lookup is a binary search.
inline std::vector<Seam> seams_as_cut(int h, int w, const std::vector<std::vector<int32_t>>& removed) {
    std::vector<std::vector<int>> cols(h, std::vector<int>(w));
    for (auto& row : cols) std::iota(row.begin(), row.end(), 0);
    std::vector<Seam> out;
    out.reserve(removed.size());
    for (const auto& ids : removed) {
        Seam s(h);
        for (int r = 0; r < h; ++r) {
            auto it = std::lower_bound(cols[r].begin(), cols[r].end(), ids[r] - r * w);
            if (it == cols[r].end() || *it != ids[r] - r * w) throw std::logic_error("seams_as_cut: pixel removed twice");
            s[r] = int(it - cols[r].begin());
            cols[r].erase(it);
        }
        out.push_back(std::move(s));
    }
    return out;
}

// Each seam as its row-0 column (little-endian uint16), then its h-1 steps four to a byte, low bits
// first: 0 = down-left, 1 = straight down, 2 = down-right. About h/4 bytes a seam, against 2*w*h bytes
// for a full int16 order map, so the visualiser can carry full-resolution photos.
inline std::vector<unsigned char> pack_seams(const std::vector<Seam>& seams) {
    std::vector<unsigned char> out;
    for (const Seam& s : seams) {
        if (s.empty() || s[0] < 0 || s[0] > 65535) throw std::invalid_argument("pack_seams: column out of range");
        out.push_back(uint8_t(s[0] & 0xff));
        out.push_back(uint8_t(s[0] >> 8));
        const size_t at = out.size();
        out.resize(at + (s.size() + 2) / 4, 0);
        for (size_t r = 1; r < s.size(); ++r) {
            int d = s[r] - s[r - 1];
            if (d < -1 || d > 1) throw std::invalid_argument("pack_seams: seam is not connected");
            out[at + (r - 1) / 4] |= uint8_t((d + 1) << (2 * ((r - 1) % 4)));
        }
    }
    return out;
}

inline std::vector<Seam> unpack_seams(const std::vector<unsigned char>& bytes, int h) {
    const size_t per = 2 + size_t(h + 2) / 4;
    if (h < 1 || bytes.size() % per) throw std::invalid_argument("unpack_seams: length is not a whole number of seams");
    std::vector<Seam> out(bytes.size() / per, Seam(h));
    for (size_t k = 0; k < out.size(); ++k) {
        const unsigned char* p = bytes.data() + per * k;
        Seam& s = out[k];
        s[0] = p[0] | p[1] << 8;
        for (int r = 1; r < h; ++r) {
            int code = (p[2 + (r - 1) / 4] >> (2 * ((r - 1) % 4))) & 3;
            if (code == 3) throw std::invalid_argument("unpack_seams: bad step code");
            s[r] = s[r - 1] + code - 1;
        }
    }
    return out;
}

// {"w":..,"h":..,"minW":..,"seams":"<base64 of pack_seams>"} for the resize slider. The browser
// replays the seams to learn which pixel each one removes; the picture itself travels separately.
inline void write_seams_json(const std::string& path, int w, int h, const std::vector<Seam>& seams) {
    ensure_parent_dir(path);
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) throw std::runtime_error("cannot write " + path);
    std::fprintf(f, "{\"w\":%d,\"h\":%d,\"minW\":%d,\"seams\":\"%s\"}\n", w, h, w - int(seams.size()),
                 base64(pack_seams(seams)).c_str());
    std::fclose(f);
}

// "300" -> 300 pixels, "50%" -> round(full * 50 / 100); the result must be in [1, full].
inline int parse_size(const std::string& s, int full) {
    size_t used = 0;
    double v = 0;
    try {
        v = std::stod(s, &used);
    } catch (const std::exception&) {
        throw std::invalid_argument("bad size '" + s + "' (expected pixels like 300 or a percentage like 50%)");
    }
    bool percent = used + 1 == s.size() && s.back() == '%';
    if (used != s.size() && !percent) throw std::invalid_argument("bad size '" + s + "'");
    int px = percent ? int(std::lround(full * v / 100.0)) : int(v);
    if (!percent && v != std::floor(v)) throw std::invalid_argument("bad size '" + s + "' (pixels must be whole)");
    if (px < 1 || px > full)
        throw std::invalid_argument("size '" + s + "' gives " + std::to_string(px) + " px, outside [1, " + std::to_string(full) + "]");
    return px;
}

}  // namespace sc

// ---------------------------------------------------------------------------------------------
// Tests (all sizes in the names are W x H, width x height; Grid/Image constructors take h, w)
// ---------------------------------------------------------------------------------------------
#ifdef SC_TESTS
namespace sc {

SC_TEST(CV1, "every method carves a 12x10 (WxH) image to exactly 7x6 with 5 + 4 seams") {
    std::mt19937 rng(601);
    Image in = random_image(rng, 10, 12);
    for (const std::string& m : finder_names()) {
        CarveOptions o;
        o.method = m;
        CarveResult r = carve(in, 7, 6, o);
        SC_CHECK_EQ(r.out.w, 7);
        SC_CHECK_EQ(r.out.h, 6);
        SC_CHECK_EQ(r.kept.w, 7);
        SC_CHECK_EQ(r.kept.h, 6);
        SC_CHECK_EQ(r.removed_v.size(), size_t(5));
        SC_CHECK_EQ(r.removed_h.size(), size_t(4));
        for (int i = 0; i < int(r.kept.a.size()); ++i) SC_CHECK_EQ(r.out.a[i], in.a[r.kept.a[i]]);
    }
}

SC_TEST(CV2, "removed_energy == sum of seam costs; first seam cost == DP on the input energy") {
    std::mt19937 rng(602);
    for (int t = 0; t < 10; ++t) {
        Image in = random_image(rng, 14, 16);
        CarveResult r = carve(in, 10, 11, CarveOptions{});
        int64_t sum = 0;
        for (int64_t c : r.seam_costs) sum += c;
        SC_CHECK_EQ(r.removed_energy, sum);
        SC_CHECK_EQ(r.seam_costs.size(), size_t(6 + 3));
        Energy e = compute_energy(to_gray(in), EnergyKind::L1);
        Stats st;
        SC_CHECK_EQ(r.seam_costs[0], seam_cost(e, dp_seam(e, st)));
    }
}

SC_TEST(CV3, "first-seam cost is identical for all 8 optimal finders (20 random 6x6)") {
    std::mt19937 rng(603);
    for (int t = 0; t < 20; ++t) {
        Image in = random_image(rng, 6, 6);
        int64_t want = -1;
        for (const std::string& m : finder_names()) {
            if (m == "greedy") continue;
            CarveOptions o;
            o.method = m;
            int64_t c = carve(in, 5, 6, o).seam_costs.at(0);
            if (want < 0) want = c;
            if (c != want) throw CheckFailure(m + ": first seam cost " + std::to_string(c) + " != " + std::to_string(want));
        }
    }
}

SC_TEST(CV4, "order-map reconstruction == real carve for every width (10 random images)") {
    std::mt19937 rng(604);
    std::uniform_int_distribution<int> dim(3, 14);
    for (int t = 0; t < 10; ++t) {
        int h = dim(rng), w = dim(rng), minW = std::max(1, w / 3);
        Image in = random_image(rng, h, w);
        IndexGrid order = removal_order(in, minW, CarveOptions{});
        for (int k = 0; k <= w - minW; ++k) SC_CHECK_GRID_EQ(width_from_order(in, order, k), carve(in, w - k, h, CarveOptions{}).out);
    }
}

SC_TEST(CV5, "overlay of a 20x15 -> 14x11 carve paints exactly 90 red and 56 blue pixels") {
    std::mt19937 rng(605);
    Image in = random_image(rng, 15, 20);
    for (RGB& p : in.a) p = {uint8_t(20 + p.r % 200), uint8_t(20 + p.g % 200), uint8_t(20 + p.b % 200)};  // no pure red/blue
    CarveResult r = carve(in, 14, 11, CarveOptions{});
    Image ov = overlay_seams(in, r);
    int red = 0, blue = 0;
    for (const RGB& p : ov.a) red += p == RGB{255, 0, 0}, blue += p == RGB{0, 0, 255};
    SC_CHECK_EQ(red, 6 * 15);
    SC_CHECK_EQ(blue, 4 * 14);
}

SC_TEST(CV6, "horizontal carve == transpose(vertical carve of transpose) for L1, L2, Gradient") {
    std::mt19937 rng(606);
    for (EnergyKind k : {EnergyKind::L1, EnergyKind::L2, EnergyKind::Gradient})
        for (int t = 0; t < 5; ++t) {
            Image in = random_image(rng, 12, 9);
            CarveOptions o;
            o.energy = k;
            CarveResult hor = carve(in, 9, 7, o);
            CarveResult ver = carve(transpose(in), 7, 9, o);
            SC_CHECK_GRID_EQ(hor.out, transpose(ver.out));
            SC_CHECK_EQ(hor.seam_costs, ver.seam_costs);
        }
}

SC_TEST(CV7, "parse_size and get_finder reject bad input") {
    SC_CHECK_EQ(parse_size("300", 800), 300);
    SC_CHECK_EQ(parse_size("50%", 801), 401);
    for (const char* bad : {"0", "900", "abc", "50%%", "12.5", "-3", "150%"}) {
        bool threw = false;
        try {
            parse_size(bad, 800);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        if (!threw) throw CheckFailure(std::string("parse_size accepted '") + bad + "'");
    }
    bool threw = false;
    try {
        get_finder("nope");
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    SC_CHECK(threw);
}

SC_TEST(CV8, "packed seams round-trip and replay to the removal order, as the slider does (10 random images)") {
    std::mt19937 rng(608);
    std::uniform_int_distribution<int> dim(1, 14);
    for (int t = 0; t < 10; ++t) {
        int h = dim(rng), w = dim(rng) + 1, minW = std::max(1, w / 3);
        Image in = random_image(rng, h, w);
        std::vector<Seam> seams = seams_as_cut(h, w, carve(in, minW, h, CarveOptions{}).removed_v);
        SC_CHECK_EQ(int(seams.size()), w - minW);
        std::vector<Seam> back = unpack_seams(pack_seams(seams), h);
        SC_CHECK_EQ(back.size(), seams.size());
        // The browser's replay: seam k picks a column of the current width, which maps to an original column.
        IndexGrid order(h, w, -1);
        std::vector<std::vector<int>> cols(h, std::vector<int>(w));
        for (auto& row : cols) std::iota(row.begin(), row.end(), 0);
        for (size_t k = 0; k < seams.size(); ++k) {
            SC_CHECK_EQ(back[k], seams[k]);
            SC_CHECK(seam_valid(seams[k], h, w - int(k)));
            for (int r = 0; r < h; ++r) {
                order.at(r, cols[r][seams[k][r]]) = int32_t(k);
                cols[r].erase(cols[r].begin() + seams[k][r]);
            }
        }
        SC_CHECK_GRID_EQ(order, removal_order(in, minW, CarveOptions{}));
    }
}

}  // namespace sc
#endif
