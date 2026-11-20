// dp.h - dynamic programming seam finder (brief 1b) and seam counting (Q2a). Owner: Vera.
// Progression, as in the lectures: exhaustive search -> memoised recursion -> bottom-up DP ->
// space and thread improvements.
//
// SRTBOT
//   Subproblem: M[i][j] = minimum total energy of a connected seam from row 0 that ends at pixel
//               (i, j) (prefix by rows, as in Unique Paths).
//   Relate:     M[i][j] = e[i][j] + min(M[i-1][j-1], M[i-1][j], M[i-1][j+1]), out of range = +inf.
//               Record the parent d in {-1, 0, +1}.
//   Topological order: increasing i (row i only needs row i-1); any column order inside a row.
//   Base:       M[0][j] = e[0][j].
//   Original:   min_j M[h-1][j], then backtrack the parents to output the seam itself.
//   Time:       h*w subproblems x O(1) relate = Theta(hw) per seam. Space Theta(hw) for the full
//               table, or Theta(w) cost rows + 1 byte per pixel of parents (rolling).
// Optimal substructure: the prefix of an optimal seam ending at (i, j) is an optimal seam ending at
// its parent, otherwise swapping in a cheaper prefix gives a cheaper seam (cut and paste).
// Overlapping subproblems: plain recursion re-solves (i, j) up to 3^(h-1-i) times; memoised and
// bottom-up DP solve each of them once.
//
// Tie rule (every variant): the leftmost parent wins on equal cost, and the leftmost end column
// wins on equal M[h-1][j]. So all DP variants return the identical seam.
#pragma once

#include <condition_variable>
#include <mutex>
#include <thread>

#include "core.h"

namespace sc {

struct DpTables {
    Grid<int64_t> M;  // h x w cumulative minimum energy
    Grid<int8_t> P;   // h x w parent offset d (row 0: 0)
    Seam seam;
    int64_t cost = 0;
};

namespace detail {

// The relate step's min with the tie rule: start with the left parent and move only when strictly
// smaller, so on equal values the leftmost parent wins.
inline int8_t pick_parent(int64_t up_left, int64_t up, int64_t up_right, int64_t& best) {
    int8_t d = -1;
    best = up_left;
    if (up < best) best = up, d = 0;
    if (up_right < best) best = up_right, d = +1;
    return d;
}

// Smallest j with minimal last[j] (the leftmost end column wins ties).
inline int argmin_first(const int64_t* last, int w) {
    return int(std::min_element(last, last + w) - last);
}

inline Seam backtrack(const Grid<int8_t>& P, int end) {
    Seam s(P.h);
    s[P.h - 1] = end;
    for (int i = P.h - 1; i > 0; --i) s[i - 1] = s[i] + P.at(i, s[i]);
    return s;
}

// Bottom-up fill. M uses the index-offset trick: width w+2 with +inf sentinels in columns 0 and
// w+1, so pixel column j is stored at j+1 and the inner loop needs no bounds checks.
inline void dp_fill(const Energy& e, Grid<int64_t>& M, Grid<int8_t>& P) {
    const int h = e.h, w = e.w;
    M = Grid<int64_t>(h, w + 2, INF);
    P = Grid<int8_t>(h, w, 0);
    for (int j = 0; j < w; ++j) M.at(0, j + 1) = e.at(0, j);
    for (int i = 1; i < h; ++i) {
        const int64_t* up = &M.at(i - 1, 0);
        int64_t* row = &M.at(i, 0);
        for (int j = 0; j < w; ++j) {
            int64_t best;
            P.at(i, j) = pick_parent(up[j], up[j + 1], up[j + 2], best);
            row[j + 1] = e.at(i, j) + best;
        }
    }
}

inline void require_nonempty(const Energy& e, const char* who) {
    if (e.h <= 0 || e.w <= 0) throw std::invalid_argument(std::string(who) + ": empty energy map");
}

}  // namespace detail

// Bottom-up DP with the full table: THE solution.
inline Seam dp_seam(const Energy& e, Stats& st) {
    detail::require_nonempty(e, "dp_seam");
    Grid<int64_t> M;
    Grid<int8_t> P;
    detail::dp_fill(e, M, P);
    st.cells += (long long)e.h * e.w;
    return detail::backtrack(P, detail::argmin_first(&M.at(e.h - 1, 1), e.w));
}

// Same as dp_seam, but returns the tables (without sentinels) for printing and tests.
inline DpTables dp_tables(const Energy& e) {
    detail::require_nonempty(e, "dp_tables");
    Grid<int64_t> Ms;
    DpTables t;
    detail::dp_fill(e, Ms, t.P);
    t.M = Grid<int64_t>(e.h, e.w);
    for (int i = 0; i < e.h; ++i)
        for (int j = 0; j < e.w; ++j) t.M.at(i, j) = Ms.at(i, j + 1);
    int end = detail::argmin_first(&Ms.at(e.h - 1, 1), e.w);
    t.seam = detail::backtrack(t.P, end);
    t.cost = t.M.at(e.h - 1, end);
    return t;
}

// Rolling rows: only two cost rows (row i-1 and row i), so Theta(w) words instead of Theta(hw).
// The 1-byte parent table is still needed to output the seam, not just its cost.
// A single row cannot be updated in place left to right: when cell j is computed, cell j-1 already
// holds row i, so j reads a "polluted value" as its up-left parent (test DP5 shows the wrong answer).
inline Seam dp_seam_rolling(const Energy& e, Stats& st) {
    detail::require_nonempty(e, "dp_seam_rolling");
    const int h = e.h, w = e.w;
    std::vector<int64_t> buf[2] = {std::vector<int64_t>(w + 2, INF), std::vector<int64_t>(w + 2, INF)};
    Grid<int8_t> P(h, w, 0);
    for (int j = 0; j < w; ++j) buf[0][j + 1] = e.at(0, j);
    for (int i = 1; i < h; ++i) {
        const std::vector<int64_t>& up = buf[(i - 1) & 1];
        std::vector<int64_t>& row = buf[i & 1];
        for (int j = 0; j < w; ++j) {
            int64_t best;
            P.at(i, j) = detail::pick_parent(up[j], up[j + 1], up[j + 2], best);
            row[j + 1] = e.at(i, j) + best;
        }
    }
    st.cells += (long long)h * w;
    return detail::backtrack(P, detail::argmin_first(&buf[(h - 1) & 1][1], w));
}

// Top-down: f(i, j) recursion with a memo table. Recursion depth = h, guarded at 2048 so it fits
// the default 1 MB Windows stack. Stats.cells counts the actual evaluations (each cell once).
inline Seam dp_seam_memo(const Energy& e, Stats& st) {
    detail::require_nonempty(e, "dp_seam_memo");
    const int h = e.h, w = e.w;
    if (h > 2048) throw std::invalid_argument("dp_seam_memo: height " + std::to_string(h) + " > 2048");
    Grid<int64_t> memo(h, w, -1);
    std::function<int64_t(int, int)> f = [&](int i, int j) -> int64_t {
        if (j < 0 || j >= w) return INF;
        int64_t& m = memo.at(i, j);
        if (m >= 0) return m;
        ++st.cells;
        if (i == 0) return m = e.at(0, j);
        int64_t best;
        detail::pick_parent(f(i - 1, j - 1), f(i - 1, j), f(i - 1, j + 1), best);
        return m = e.at(i, j) + best;
    };
    std::vector<int64_t> last(w);
    for (int j = 0; j < w; ++j) last[j] = f(h - 1, j);
    // Backtrack by recomputing the argmin over the memoised row above.
    Seam s(h);
    s[h - 1] = detail::argmin_first(last.data(), w);
    for (int i = h - 1; i > 0; --i) {
        int j = s[i];
        int64_t best;
        s[i - 1] = j + detail::pick_parent(f(i - 1, j - 1), f(i - 1, j), f(i - 1, j + 1), best);
    }
    return s;
}

namespace detail {

// Reusable barrier (std::barrier is C++20): the last thread to arrive starts a new generation.
class Barrier {
public:
    explicit Barrier(int n) : n_(n) {}
    void wait() {
        std::unique_lock<std::mutex> lock(m_);
        long long gen = gen_;
        if (++count_ == n_) {
            count_ = 0;
            ++gen_;
            cv_.notify_all();
        } else {
            cv_.wait(lock, [&] { return gen_ != gen; });
        }
    }

private:
    std::mutex m_;
    std::condition_variable cv_;
    int n_, count_ = 0;
    long long gen_ = 0;
};

}  // namespace detail

// Row by row, with the columns of a row split across threads. Workers start once per call, each
// owns a fixed column chunk, and all meet at a barrier after every row (row i needs all of row i-1).
// That is h barriers for only ~w operations per row, so expect it to be slower than dp_seam at
// these sizes: parallelise the energy map, not the seam DP. Same tie rule, so the identical seam.
inline Seam dp_seam_mt(const Energy& e, Stats& st, int threads) {
    detail::require_nonempty(e, "dp_seam_mt");
    const int h = e.h, w = e.w;
    const int t = std::clamp(threads, 1, w);
    Grid<int64_t> M(h, w + 2, INF);
    Grid<int8_t> P(h, w, 0);
    for (int j = 0; j < w; ++j) M.at(0, j + 1) = e.at(0, j);
    detail::Barrier barrier(t);
    auto worker = [&](int k) {
        const int j0 = int(int64_t(w) * k / t), j1 = int(int64_t(w) * (k + 1) / t);
        for (int i = 1; i < h; ++i) {
            const int64_t* up = &M.at(i - 1, 0);
            int64_t* row = &M.at(i, 0);
            for (int j = j0; j < j1; ++j) {
                int64_t best;
                P.at(i, j) = detail::pick_parent(up[j], up[j + 1], up[j + 2], best);
                row[j + 1] = e.at(i, j) + best;
            }
            barrier.wait();
        }
    };
    std::vector<std::thread> pool;
    for (int k = 1; k < t; ++k) pool.emplace_back(worker, k);
    worker(0);
    for (auto& th : pool) th.join();
    st.cells += (long long)h * w;
    return detail::backtrack(P, detail::argmin_first(&M.at(h - 1, 1), w));
}

// Exhaustive search over all seams (w * 3^(h-1) at most): the correctness oracle.
// Seams are enumerated in lexicographic order and only a strictly cheaper seam replaces the best,
// so ties go to the lexicographically first seam. Stats.cells = number of complete seams examined.
inline Seam brute_seam(const Energy& e, Stats& st) {
    detail::require_nonempty(e, "brute_seam");
    const int h = e.h, w = e.w;
    if (double(w) * std::pow(3.0, h - 1) > 5e7)
        throw std::invalid_argument("brute_seam: too many seams for " + std::to_string(h) + "x" + std::to_string(w));
    Seam cur(h), best;
    int64_t best_cost = INF;
    std::function<void(int, int64_t)> go = [&](int i, int64_t cost) {
        if (i == h) {
            ++st.cells;
            if (cost < best_cost) best_cost = cost, best = cur;
            return;
        }
        int lo = i == 0 ? 0 : std::max(cur[i - 1] - 1, 0);
        int hi = i == 0 ? w - 1 : std::min(cur[i - 1] + 1, w - 1);
        for (int j = lo; j <= hi; ++j) {
            cur[i] = j;
            go(i + 1, cost + e.at(i, j));
        }
    };
    go(0, 0);
    return best;
}

// Q2a(i): number of vertical seams. N[0][j] = 1, N[i][j] = sum over d of N[i-1][j+d].
// Saturates at UINT64_MAX instead of overflowing.
inline uint64_t count_seams(int h, int w) {
    if (h <= 0 || w <= 0) return 0;
    auto add = [](uint64_t a, uint64_t b) { return a > UINT64_MAX - b ? UINT64_MAX : a + b; };
    std::vector<uint64_t> prev(w, 1), cur(w);
    for (int i = 1; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            uint64_t n = prev[j];
            if (j > 0) n = add(n, prev[j - 1]);
            if (j + 1 < w) n = add(n, prev[j + 1]);
            cur[j] = n;
        }
        prev.swap(cur);
    }
    uint64_t total = 0;
    for (uint64_t n : prev) total = add(total, n);
    return total;
}

// Enumerates every seam (small h and w only): the test oracle for count_seams.
inline uint64_t count_seams_brute(int h, int w) {
    std::function<uint64_t(int, int)> from = [&](int i, int j) -> uint64_t {
        if (j < 0 || j >= w) return 0;
        if (i == h - 1) return 1;
        return from(i + 1, j - 1) + from(i + 1, j) + from(i + 1, j + 1);
    };
    uint64_t total = 0;
    for (int j = 0; j < w; ++j) total += from(0, j);
    return total;
}

// Pretty-prints E, M and P (parent arrows: \ up-left, | up, / up-right) with the seam cells
// marked '*', then the seam and its cost. Used by `seamcarve dumpdp` for the presentation.
inline std::string format_dp_tables(const Energy& e, const DpTables& t) {
    std::ostringstream o;
    char buf[32];
    auto header = [&](const char* title) {
        o << title << "\n     |";
        for (int j = 0; j < e.w; ++j) std::snprintf(buf, sizeof buf, j ? " %5d" : "%5d", j), o << buf;
        o << "\n";
    };
    auto rows = [&](auto cell) {
        for (int i = 0; i < e.h; ++i) {
            std::string line;
            std::snprintf(buf, sizeof buf, "%4d |", i);
            line += buf;
            for (int j = 0; j < e.w; ++j) line += cell(i, j) + (t.seam[i] == j ? "*" : " ");
            while (!line.empty() && line.back() == ' ') line.pop_back();
            o << line << "\n";
        }
        o << "\n";
    };
    auto num = [&](long long v) {
        std::snprintf(buf, sizeof buf, "%5lld", v);
        return std::string(buf);
    };
    header("E (energy)");
    rows([&](int i, int j) { return num(e.at(i, j)); });
    header("M (minimum seam cost ending here)");
    rows([&](int i, int j) { return num(t.M.at(i, j)); });
    header("P (parent: \\ up-left, | up, / up-right)");
    rows([&](int i, int j) {
        const char* arrow = i == 0 ? "." : t.P.at(i, j) < 0 ? "\\" : t.P.at(i, j) == 0 ? "|" : "/";
        return std::string("    ") + arrow;
    });
    o << "seam " << t.seam << "  cost " << t.cost << "\n";
    return o.str();
}

inline Energy energy_from(const std::vector<std::vector<int>>& rows) {
    Energy e(int(rows.size()), int(rows[0].size()));
    for (int i = 0; i < e.h; ++i)
        for (int j = 0; j < e.w; ++j) e.at(i, j) = rows[i][j];
    return e;
}

// The hand-checked 4x5 example (test DP6), printed by `seamcarve dumpdp --demo`.
inline Energy dp_demo_energy() {
    return energy_from({{3, 1, 4, 1, 5}, {9, 2, 6, 5, 3}, {5, 8, 9, 7, 9}, {3, 2, 3, 8, 4}});
}

}  // namespace sc

// ---------------------------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------------------------
#ifdef SC_TESTS
namespace sc {

SC_TEST(DP1, "DP cost == brute-force cost (2,000 random grids up to 8x8, many ties)") {
    std::mt19937 rng(301);
    std::uniform_int_distribution<int> dim(1, 8);
    for (int t = 0; t < 2000; ++t) {
        Energy e = random_energy(rng, dim(rng), dim(rng), t % 2 ? 3 : 2040);
        Stats st;
        SC_CHECK_EQ(seam_cost(e, dp_seam(e, st)), seam_cost(e, brute_seam(e, st)));
    }
}

SC_TEST(DP2, "DP seam is valid and its cost == min of the last row of M") {
    std::mt19937 rng(302);
    std::uniform_int_distribution<int> dim(1, 30);
    for (int t = 0; t < 300; ++t) {
        Energy e = random_energy(rng, dim(rng), dim(rng), 2040);
        Stats st;
        Seam s = dp_seam(e, st);
        std::string why;
        SC_CHECK(seam_valid(s, e.h, e.w, &why));
        DpTables tb = dp_tables(e);
        int64_t last_min = *std::min_element(tb.M.a.end() - e.w, tb.M.a.end());
        SC_CHECK_EQ(seam_cost(e, s), last_min);
        SC_CHECK_EQ(tb.cost, last_min);
        SC_CHECK_EQ(tb.seam, s);
    }
}

SC_TEST(DP3, "memoised recursion gives the identical seam to bottom-up (500 random up to 40x40)") {
    std::mt19937 rng(303);
    std::uniform_int_distribution<int> dim(1, 40);
    for (int t = 0; t < 500; ++t) {
        Energy e = random_energy(rng, dim(rng), dim(rng), t % 2 ? 3 : 2040);
        Stats a, b;
        SC_CHECK_EQ(dp_seam_memo(e, a), dp_seam(e, b));
        SC_CHECK_EQ(a.cells, b.cells);  // each subproblem is evaluated exactly once
    }
}

SC_TEST(DP4, "rolling and 3-thread DP give the identical seam (500 random up to 60x60)") {
    std::mt19937 rng(304);
    std::uniform_int_distribution<int> dim(1, 60);
    for (int t = 0; t < 500; ++t) {
        Energy e = random_energy(rng, dim(rng), t < 20 ? 1 + t % 3 : dim(rng), t % 2 ? 3 : 2040);
        Stats st;
        Seam want = dp_seam(e, st);
        SC_CHECK_EQ(dp_seam_rolling(e, st), want);
        SC_CHECK_EQ(dp_seam_mt(e, st, 3), want);
    }
}

SC_TEST(DP5, "one-row in-place DP reads a polluted value: 9 instead of the true 6") {
    Energy e = energy_from({{5, 8, 9}, {9, 1, 3}});
    // Deliberately wrong: update a single row left to right, so cell j-1 already holds row 1.
    std::vector<int64_t> row(e.a.begin(), e.a.begin() + e.w);
    for (int i = 1; i < e.h; ++i)
        for (int j = 0; j < e.w; ++j) {
            int64_t best = row[j];
            if (j > 0) best = std::min(best, row[j - 1]);
            if (j + 1 < e.w) best = std::min(best, row[j + 1]);
            row[j] = e.at(i, j) + best;
        }
    SC_CHECK_EQ(*std::min_element(row.begin(), row.end()), int64_t(9));
    Stats st;
    SC_CHECK_EQ(seam_cost(e, dp_seam(e, st)), int64_t(6));
}

SC_TEST(DP6, "hand-filled 4x5 table: M, P, seam [1,1,0,1], cost 10") {
    DpTables t = dp_tables(dp_demo_energy());
    Grid<int64_t> M(4, 5);
    M.a = {3, 1, 4, 1, 5, 10, 3, 7, 6, 4, 8, 11, 12, 11, 13, 11, 10, 14, 19, 15};
    Grid<int8_t> P(4, 5);
    P.a = {0, 0, 0, 0, 0, 1, 0, -1, 0, -1, 1, 0, -1, 1, 0, 0, -1, -1, 0, -1};
    SC_CHECK_GRID_EQ(t.M, M);
    SC_CHECK_GRID_EQ(t.P, P);
    SC_CHECK_EQ(t.seam, (Seam{1, 1, 0, 1}));
    SC_CHECK_EQ(t.cost, int64_t(10));
    Stats st;
    SC_CHECK_EQ(dp_seam(dp_demo_energy(), st), (Seam{1, 1, 0, 1}));
}

SC_TEST(DP7, "edge cases: w=1, h=1 (leftmost tie), all-equal grid") {
    Stats st;
    SC_CHECK_EQ(dp_seam(Energy(5, 1, 7), st), Seam(5, 0));
    SC_CHECK_EQ(dp_seam(energy_from({{4, 2, 9, 2}}), st), Seam{1});
    SC_CHECK_EQ(dp_seam(Energy(6, 6, 3), st), Seam(6, 0));
    SC_CHECK_EQ(dp_seam_memo(Energy(6, 6, 3), st), Seam(6, 0));
    SC_CHECK_EQ(dp_seam_rolling(Energy(1, 1, 0), st), Seam{0});
    SC_CHECK_EQ(brute_seam(Energy(1, 1, 0), st), Seam{0});
}

SC_TEST(DP8, "count_seams == enumeration; (3,3) = 17; at least 2^(h-1) for w >= 2") {
    for (int h = 1; h <= 7; ++h)
        for (int w = 1; w <= 6; ++w) SC_CHECK_EQ(count_seams(h, w), count_seams_brute(h, w));
    SC_CHECK_EQ(count_seams(3, 3), uint64_t(17));
    for (int h = 1; h <= 30; ++h)
        for (int w : {2, 3, 10, 100}) SC_CHECK(count_seams(h, w) >= (uint64_t(1) << (h - 1)));
    SC_CHECK_EQ(count_seams(200, 200), UINT64_MAX);  // saturates instead of overflowing
}

SC_TEST(DP9, "format_dp_tables shows the DP6 row 3 of M and the seam") {
    Energy e = dp_demo_energy();
    std::string out = format_dp_tables(e, dp_tables(e));
    SC_CHECK(out.find("   3 |   11    10*   14    19    15\n") != std::string::npos);
    SC_CHECK(out.find("seam [1,1,0,1]  cost 10") != std::string::npos);
}

}  // namespace sc
#endif
