// dp.h - dynamic programming seam finder (brief 1b). Owner: Vera.
// Exhaustive search is the correctness oracle; bottom-up DP is the solution.
//
// SRTBOT
//   Subproblem: M[i][j] = minimum total energy of a connected seam from row 0 that ends at pixel
//               (i, j) (prefix by rows, as in Unique Paths).
//   Relate:     M[i][j] = e[i][j] + min(M[i-1][j-1], M[i-1][j], M[i-1][j+1]), out of range = +inf.
//               Record the parent d in {-1, 0, +1}.
//   Topological order: increasing i (row i only needs row i-1); any column order inside a row.
//   Base:       M[0][j] = e[0][j].
//   Original:   min_j M[h-1][j], then backtrack the parents to output the seam itself.
//   Time:       h*w subproblems x O(1) relate = Theta(hw) per seam. Space Theta(hw) for the table.
// Optimal substructure: the prefix of an optimal seam ending at (i, j) is an optimal seam ending at
// its parent, otherwise swapping in a cheaper prefix gives a cheaper seam (cut and paste).
// Overlapping subproblems: plain recursion re-solves (i, j) up to 3^(h-1-i) times; bottom-up DP
// solves each of them once.
//
// Tie rule: the leftmost parent wins on equal cost, and the leftmost end column wins on equal
// M[h-1][j].
#pragma once

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

inline Energy energy_from(const std::vector<std::vector<int>>& rows) {
    Energy e(int(rows.size()), int(rows[0].size()));
    for (int i = 0; i < e.h; ++i)
        for (int j = 0; j < e.w; ++j) e.at(i, j) = rows[i][j];
    return e;
}

// The hand-checked 4x5 example (test DP6).
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
    SC_CHECK_EQ(brute_seam(Energy(1, 1, 0), st), Seam{0});
}

}  // namespace sc
#endif
