// graphcut.h - seam finding as a minimum s-t cut, solved by maximum flow (brief 4). Owner: Zafir.
// Solvers from the max-flow lecture: Ford-Fulkerson (DFS), Edmonds-Karp (BFS) and Dinic.
//
// Graph model for a vertical seam on an h x w energy map:
//   - Pad one dummy column (energy 0) on the right, so W = w + 1. Without it the last real column
//     can never be cut (test GC3: unpadded min cut != optimum on [[9,9,0],[9,9,0]]).
//   - Nodes: pixel (i, j) -> i*W + j; s = h*W, t = h*W + 1.
//   - Edges:
//       s -> (i, 0)                INF
//       (i, W-1) -> t              INF
//       (i, j) -> (i, j+1)         e[i][j]   cutting it = deleting pixel (i, j)       (j < W-1)
//       (i, j+1) -> (i, j)         INF       makes each row's source side a prefix    (j < W-1)
//       (i, j) -> (i-1, j-1)       INF       if i > 0      } force |seam[i] - seam[i+-1]| <= 1
//       (i, j) -> (i+1, j-1)       INF       if i < h-1    }                           (j >= 1)
//     Edge pairs: h(2 + 2(W-1)) + 2(h-1)(W-1)  (test GC1).
//   - Seam: seam[i] = the rightmost j with (i, j) in S, the source side of the minimum cut.
//     Column 0 is always in S; the dummy column never is (its INF edge to t would cut at INF).
// Correctness (max-flow min-cut theorem): every finite cut crosses each row's forward chain exactly
// once, at pixels forming a connected seam, so cut capacity = seam cost. Minimum cut = optimal seam
// = maximum flow value.
// Which edges are essential? Let k_i be the rightmost S-pixel of row i. The diagonal INF edges
// alone force |k_i - k_{i+-1}| <= 1, and any finite cut pays at least e[i][k_i] per row, so even
// without the backward INF edges the min-cut value and the rightmost-in-S seam stay optimal. The
// backward edges make S a per-row prefix, so every finite cut IS exactly one seam (a clean
// one-to-one correspondence). Only the structural tests GC1/GC2 notice if they are missing.
//
// Complexity: n = hw + 2 nodes, m ~ 5hw edge pairs. FF is O(f m) with f up to 2040h,
// Edmonds-Karp O(n m^2), Dinic O(n^2 m): all far slower than DP's Theta(hw) for the same seam.
#pragma once

#include <tuple>

#include "core.h"

namespace sc {

// Adjacency list stored as edge arrays. Edge e and its reverse e^1 are added together.
struct FlowGraph {
    int n = 0, s = 0, t = 0;
    std::vector<int> head, to, nxt;   // head[u] = first edge of u (-1: none), nxt[e] = next edge of the same tail
    std::vector<int64_t> cap, cap0;   // residual and original capacity

    explicit FlowGraph(int n_ = 0) : n(n_), head(n_, -1) {}

    int add_edge(int u, int v, int64_t c) {
        int e = int(to.size());
        for (auto [a, b, k] : {std::tuple{u, v, c}, std::tuple{v, u, int64_t(0)}}) {
            to.push_back(b);
            nxt.push_back(head[a]);
            head[a] = int(to.size()) - 1;
            cap.push_back(k);
            cap0.push_back(k);
        }
        return e;
    }
    int tail(int e) const { return to[e ^ 1]; }
};

struct SeamGraph {
    FlowGraph g;
    int h = 0, w = 0, W = 0;  // W = w + 1 when padded
    int node(int i, int j) const { return i * W + j; }
};

inline SeamGraph build_seam_graph(const Energy& e, bool pad = true) {
    if (e.h <= 0 || e.w <= 0) throw std::invalid_argument("build_seam_graph: empty energy map");
    SeamGraph sg;
    sg.h = e.h;
    sg.w = e.w;
    sg.W = pad ? e.w + 1 : e.w;
    FlowGraph& g = sg.g;
    g = FlowGraph(sg.h * sg.W + 2);
    g.s = sg.h * sg.W;
    g.t = g.s + 1;
    for (int i = 0; i < sg.h; ++i) {
        g.add_edge(g.s, sg.node(i, 0), INF);
        g.add_edge(sg.node(i, sg.W - 1), g.t, INF);
        for (int j = 0; j + 1 < sg.W; ++j) {
            g.add_edge(sg.node(i, j), sg.node(i, j + 1), e.at(i, j));  // j < w here, also when padded
            g.add_edge(sg.node(i, j + 1), sg.node(i, j), INF);
        }
        for (int j = 1; j < sg.W; ++j) {
            if (i > 0) g.add_edge(sg.node(i, j), sg.node(i - 1, j - 1), INF);
            if (i + 1 < sg.h) g.add_edge(sg.node(i, j), sg.node(i + 1, j - 1), INF);
        }
    }
    return sg;
}

namespace detail {

// Pushes x units along the path: the forward residual drops and the reverse residual grows, which
// is what lets a later path "undo" flow (the key idea of Ford-Fulkerson).
inline void augment(FlowGraph& g, const std::vector<int>& path, int64_t x) {
    for (int e : path) {
        g.cap[e] -= x;
        g.cap[e ^ 1] += x;
    }
}

inline int64_t bottleneck(const FlowGraph& g, const std::vector<int>& path) {
    int64_t x = INF;
    for (int e : path) x = std::min(x, g.cap[e]);
    if (x >= INF) throw std::logic_error("augmenting path of infinite capacity");
    return x;
}

// Path from s to v as edge ids, following parent edges back from v.
inline std::vector<int> path_to(const FlowGraph& g, const std::vector<int>& parent, int v) {
    std::vector<int> path;
    for (; v != g.s; v = g.tail(parent[v])) path.push_back(parent[v]);
    std::reverse(path.begin(), path.end());
    return path;
}

}  // namespace detail

// Ford-Fulkerson: any augmenting path, found by DFS with an explicit stack (paths can be h*W long,
// too deep for recursion). O(f m).
inline int64_t maxflow_ff(FlowGraph& g, Stats& st) {
    int64_t flow = 0;
    std::vector<int> parent(g.n), stack;
    std::vector<char> seen(g.n);
    while (true) {
        std::fill(seen.begin(), seen.end(), 0);
        stack.assign(1, g.s);
        seen[g.s] = 1;
        while (!stack.empty() && !seen[g.t]) {
            int u = stack.back();
            stack.pop_back();
            ++st.bfs_visits;
            for (int e = g.head[u]; e != -1; e = g.nxt[e])
                if (g.cap[e] > 0 && !seen[g.to[e]]) {
                    seen[g.to[e]] = 1;
                    parent[g.to[e]] = e;
                    stack.push_back(g.to[e]);
                }
        }
        if (!seen[g.t]) return flow;
        std::vector<int> path = detail::path_to(g, parent, g.t);
        int64_t x = detail::bottleneck(g, path);
        detail::augment(g, path, x);
        flow += x;
        ++st.augment;
    }
}

// Edmonds-Karp: always the shortest augmenting path (BFS). O(n m^2).
inline int64_t maxflow_ek(FlowGraph& g, Stats& st) {
    int64_t flow = 0;
    std::vector<int> parent(g.n), queue;
    std::vector<char> seen(g.n);
    while (true) {
        std::fill(seen.begin(), seen.end(), 0);
        queue.assign(1, g.s);
        seen[g.s] = 1;
        for (size_t q = 0; q < queue.size() && !seen[g.t]; ++q) {
            int u = queue[q];
            ++st.bfs_visits;
            for (int e = g.head[u]; e != -1; e = g.nxt[e])
                if (g.cap[e] > 0 && !seen[g.to[e]]) {
                    seen[g.to[e]] = 1;
                    parent[g.to[e]] = e;
                    queue.push_back(g.to[e]);
                }
        }
        if (!seen[g.t]) return flow;
        std::vector<int> path = detail::path_to(g, parent, g.t);
        int64_t x = detail::bottleneck(g, path);
        detail::augment(g, path, x);
        flow += x;
        ++st.augment;
    }
}

// Dinic: per phase, BFS levels from s, then a blocking flow using only edges that go one level
// deeper. The DFS is iterative (paths can be long) and keeps a current-arc pointer it[u], so each
// edge is skipped at most once per phase. O(n^2 m).
inline int64_t maxflow_dinic(FlowGraph& g, Stats& st) {
    int64_t flow = 0;
    std::vector<int> level(g.n), it(g.n), queue, path;
    while (true) {
        std::fill(level.begin(), level.end(), -1);
        level[g.s] = 0;
        queue.assign(1, g.s);
        for (size_t q = 0; q < queue.size(); ++q) {
            int u = queue[q];
            ++st.bfs_visits;
            for (int e = g.head[u]; e != -1; e = g.nxt[e])
                if (g.cap[e] > 0 && level[g.to[e]] < 0) {
                    level[g.to[e]] = level[u] + 1;
                    queue.push_back(g.to[e]);
                }
        }
        if (level[g.t] < 0) return flow;
        ++st.phases;
        it = g.head;
        path.clear();
        int u = g.s;
        while (true) {
            if (u == g.t) {
                int64_t x = detail::bottleneck(g, path);
                detail::augment(g, path, x);
                flow += x;
                ++st.augment;
                // Retreat to the tail of the first saturated edge and continue from there.
                size_t k = 0;
                while (g.cap[path[k]] > 0) ++k;
                u = g.tail(path[k]);
                path.resize(k);
                continue;
            }
            int& e = it[u];
            while (e != -1 && !(g.cap[e] > 0 && level[g.to[e]] == level[u] + 1)) e = g.nxt[e];
            if (e != -1) {
                path.push_back(e);
                u = g.to[e];
            } else {
                if (u == g.s) break;  // blocking flow reached
                level[u] = -1;        // dead end: never enter u again in this phase
                int back = path.back();
                path.pop_back();
                u = g.tail(back);
                it[u] = g.nxt[it[u]];
            }
        }
    }
}

// Nodes reachable from s in the residual graph: the source side S of a minimum cut.
inline std::vector<char> source_side(const FlowGraph& g) {
    std::vector<char> S(g.n, 0);
    std::vector<int> queue(1, g.s);
    S[g.s] = 1;
    for (size_t q = 0; q < queue.size(); ++q)
        for (int e = g.head[queue[q]]; e != -1; e = g.nxt[e])
            if (g.cap[e] > 0 && !S[g.to[e]]) {
                S[g.to[e]] = 1;
                queue.push_back(g.to[e]);
            }
    return S;
}

inline Seam seam_from_cut(const SeamGraph& sg, const std::vector<char>& S) {
    Seam s(sg.h, -1);
    for (int i = 0; i < sg.h; ++i)
        for (int j = 0; j < sg.W; ++j)
            if (S[sg.node(i, j)]) s[i] = j;  // rightmost pixel of row i in S
    return s;
}

struct FlowCheck {
    bool ok = false;
    std::string why;
    int64_t value = 0;
};

// Certificate that the residual graph holds a maximum flow, using only cap0 and cap (so it does
// not trust the solver): capacity and skew symmetry per edge, conservation, maximality (t not
// reachable) and min cut == flow value.
inline FlowCheck verify_flow(const FlowGraph& g) {
    FlowCheck r;
    auto fail = [&](const std::string& why) {
        r.ok = false;
        r.why = why;
        return r;
    };
    std::vector<int64_t> net(g.n, 0);  // flow out minus flow in
    for (int e = 0; e < int(g.to.size()); e += 2) {
        int64_t f = g.cap0[e] - g.cap[e];
        if (f < 0 || f > g.cap0[e]) return fail("edge " + std::to_string(e) + ": flow " + std::to_string(f) + " outside capacity");
        if (g.cap[e ^ 1] != g.cap0[e ^ 1] + f)
            return fail("edge " + std::to_string(e) + ": reverse residual " + std::to_string(g.cap[e ^ 1]) +
                        " != flow " + std::to_string(f));
        net[g.tail(e)] += f;
        net[g.to[e]] -= f;
    }
    for (int v = 0; v < g.n; ++v)
        if (v != g.s && v != g.t && net[v] != 0) return fail("node " + std::to_string(v) + ": flow not conserved");
    r.value = net[g.s];
    std::vector<char> S = source_side(g);
    if (S[g.t]) return fail("t is reachable in the residual graph: the flow is not maximum");
    int64_t cut = 0;
    for (int e = 0; e < int(g.to.size()); e += 2)
        if (S[g.tail(e)] && !S[g.to[e]]) cut += g.cap0[e];
    if (cut != r.value) return fail("cut " + std::to_string(cut) + " != flow " + std::to_string(r.value));
    r.ok = true;
    return r;
}

namespace detail {

template <class Solver>
Seam cut_seam(const Energy& e, Stats& st, Solver solve) {
    SeamGraph sg = build_seam_graph(e);
    solve(sg.g, st);
    st.cells += sg.g.n;
    return seam_from_cut(sg, source_side(sg.g));
}

}  // namespace detail

inline Seam ff_seam(const Energy& e, Stats& st) { return detail::cut_seam(e, st, maxflow_ff); }
inline Seam ek_seam(const Energy& e, Stats& st) { return detail::cut_seam(e, st, maxflow_ek); }
inline Seam dinic_seam(const Energy& e, Stats& st) { return detail::cut_seam(e, st, maxflow_dinic); }

}  // namespace sc

// ---------------------------------------------------------------------------------------------
// Tests (brute force and DP are the oracles)
// ---------------------------------------------------------------------------------------------
#ifdef SC_TESTS
#include "dp.h"

namespace sc {

using MaxFlow = int64_t (*)(FlowGraph&, Stats&);
const std::pair<const char*, MaxFlow> kSolvers[] = {{"ff", maxflow_ff}, {"ek", maxflow_ek}, {"dinic", maxflow_dinic}};

SC_TEST(GC1, "edge-pair count == h(2 + 2(W-1)) + 2(h-1)(W-1), padded and unpadded") {
    for (auto [h, w] : {std::pair{1, 1}, std::pair{1, 5}, std::pair{3, 1}, std::pair{4, 6}, std::pair{7, 7}})
        for (bool pad : {true, false}) {
            SeamGraph sg = build_seam_graph(Energy(h, w, 1), pad);
            int W = pad ? w + 1 : w;
            SC_CHECK_EQ(sg.W, W);
            SC_CHECK_EQ(int(sg.g.to.size()) / 2, h * (2 + 2 * (W - 1)) + 2 * (h - 1) * (W - 1));
        }
}

SC_TEST(GC2, "padded 2x2 graph has exactly the listed edges and capacities") {
    Energy e(2, 2);
    e.a = {10, 11, 20, 21};
    SeamGraph sg = build_seam_graph(e);
    // W = 3: (0,0)=0 (0,1)=1 (0,2)=2 (1,0)=3 (1,1)=4 (1,2)=5, s=6, t=7.
    std::vector<std::tuple<int, int, int64_t>> want = {
        {6, 0, INF}, {6, 3, INF}, {2, 7, INF}, {5, 7, INF},                    // source and sink
        {0, 1, 10}, {1, 2, 11}, {3, 4, 20}, {4, 5, 21},                        // forward chains
        {1, 0, INF}, {2, 1, INF}, {4, 3, INF}, {5, 4, INF},                    // backward INF
        {4, 0, INF}, {5, 1, INF}, {1, 3, INF}, {2, 4, INF}};                   // diagonals
    std::vector<std::tuple<int, int, int64_t>> got;
    const FlowGraph& g = sg.g;
    for (int k = 0; k < int(g.to.size()); k += 2) {
        got.emplace_back(g.tail(k), g.to[k], g.cap0[k]);
        SC_CHECK_EQ(g.tail(k ^ 1), g.to[k]);  // the reverse edge goes back, with capacity 0
        SC_CHECK_EQ(g.cap0[k ^ 1], int64_t(0));
    }
    std::sort(want.begin(), want.end());
    std::sort(got.begin(), got.end());
    SC_CHECK(got == want);
}

SC_TEST(GC3, "padding is needed: [[9,9,0],[9,9,0]] has optimum 0, unpadded flow != 0, padded == 0") {
    Energy e(2, 3);
    e.a = {9, 9, 0, 9, 9, 0};
    Stats st;
    SC_CHECK_EQ(seam_cost(e, dp_seam(e, st)), int64_t(0));
    for (auto [name, solve] : kSolvers) {
        SeamGraph unpadded = build_seam_graph(e, false), padded = build_seam_graph(e, true);
        SC_CHECK(solve(unpadded.g, st) != 0);
        SC_CHECK_EQ(solve(padded.g, st), int64_t(0));
    }
}

// The GC4/GC5 grids: 500 random, h, w in 1..7, half of them with values 0..3 (many ties).
template <class F>
void for_flow_grids(F f) {
    std::mt19937 rng(504);
    std::uniform_int_distribution<int> dim(1, 7);
    for (int t = 0; t < 500; ++t) f(random_energy(rng, dim(rng), dim(rng), t % 2 ? 3 : 2040));
}

SC_TEST(GC4, "FF == EK == Dinic == brute-force optimum (500 random grids)") {
    for_flow_grids([](const Energy& e) {
        Stats st;
        int64_t best = seam_cost(e, brute_seam(e, st));
        for (auto [name, solve] : kSolvers) {
            SeamGraph sg = build_seam_graph(e);
            SC_CHECK_EQ(solve(sg.g, st), best);
        }
    });
}

SC_TEST(GC5, "flow certificate holds for every solver, and fails on a tampered flow") {
    int tampered_rejected = 0;
    for_flow_grids([&](const Energy& e) {
        for (auto [name, solve] : kSolvers) {
            SeamGraph sg = build_seam_graph(e);
            Stats st;
            int64_t v = solve(sg.g, st);
            FlowCheck c = verify_flow(sg.g);
            if (!c.ok) throw CheckFailure(std::string(name) + ": " + c.why);
            SC_CHECK_EQ(c.value, v);
            // Tamper with one forward chain edge (edge pair 2 is (0,0) -> (0,1)): +1 unit of "flow".
            sg.g.cap[4] -= 1;
            tampered_rejected += !verify_flow(sg.g).ok;
        }
    });
    SC_CHECK_EQ(tampered_rejected, 3 * 500);
}

SC_TEST(GC6, "the cut seam is valid and its cost == the flow value (every solver)") {
    std::mt19937 rng(506);
    std::uniform_int_distribution<int> dim(1, 12);
    for (int t = 0; t < 300; ++t) {
        Energy e = random_energy(rng, dim(rng), dim(rng), t % 2 ? 3 : 2040);
        for (auto [name, solve] : kSolvers) {
            SeamGraph sg = build_seam_graph(e);
            Stats st;
            int64_t v = solve(sg.g, st);
            Seam s = seam_from_cut(sg, source_side(sg.g));
            std::string why;
            if (!seam_valid(s, e.h, e.w, &why)) throw CheckFailure(std::string(name) + ": " + why);
            SC_CHECK_EQ(seam_cost(e, s), v);
        }
    }
}

SC_TEST(GC7, "32x32: Edmonds-Karp and Dinic seams cost the same as DP") {
    std::mt19937 rng(507);
    for (int t = 0; t < 3; ++t) {
        Energy e = random_energy(rng, 32, 32, 2040);
        Stats st;
        int64_t d = seam_cost(e, dp_seam(e, st));
        SC_CHECK_EQ(seam_cost(e, dinic_seam(e, st)), d);
        SC_CHECK_EQ(seam_cost(e, ek_seam(e, st)), d);
    }
}

SC_TEST(GC8, "edge cases: w=1 gives all zeros; h=1 gives the minimum column") {
    std::mt19937 rng(508);
    for (auto solve : {ff_seam, ek_seam, dinic_seam}) {
        Stats st;
        SC_CHECK_EQ(solve(random_energy(rng, 6, 1, 2040), st), Seam(6, 0));
        Energy row = random_energy(rng, 1, 9, 2040);
        SC_CHECK_EQ(seam_cost(row, solve(row, st)), int64_t(*std::min_element(row.a.begin(), row.a.end())));
        SC_CHECK_EQ(solve(Energy(1, 1, 5), st), Seam{0});
    }
}

SC_TEST(GC9, "work counters are filled; FF needs at least as many augmenting paths as EK somewhere") {
    std::mt19937 rng(509);
    bool ff_at_least_ek = false;
    for (int t = 0; t < 5; ++t) {
        Energy e = random_energy(rng, 16, 16, 2040);
        Stats ff, ek, dinic;
        ff_seam(e, ff);
        ek_seam(e, ek);
        dinic_seam(e, dinic);
        for (const Stats* s : {&ff, &ek, &dinic}) {
            SC_CHECK(s->augment > 0);
            SC_CHECK(s->bfs_visits > 0);
            SC_CHECK(s->cells > 0);
        }
        SC_CHECK(dinic.phases > 0);
        ff_at_least_ek |= ff.augment >= ek.augment;
    }
    SC_CHECK(ff_at_least_ek);
}

}  // namespace sc
#endif
