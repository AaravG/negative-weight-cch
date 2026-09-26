// Potential-free Customizable Contraction Hierarchies for negative edge weights.
//
//   cch <graph-dir> [options]
//     --order <file>     read/write the nested-dissection order (binary int32)
//     --queries <n>      number of random queries (default 1000)
//     --check <n>        verify this many of them against Dijkstra + Johnson potential
//     --metric <name>    metric to use (default: all in meta.txt)
//     --shifted-cch      also run the classical pipelines on the same hierarchy: shift by a
//                        height-induced or Johnson potential, customize the shifted metric,
//                        query it with the elimination-tree query and with the Dijkstra-based
//                        CCH query (with and without stall-on-demand)
//     --parallel <t>     also run level-parallel customization with t threads
//     --repeat <r>       repetitions per timing; medians are reported (default 5)
//
// Graph directory: meta.txt (n, m, metric names), off.i64, tgt.i32, x.f64, y.f64,
// w_<metric>.f64 - written by export_graph.py.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <deque>
#include <filesystem>
#include <fstream>
#include <limits>
#include <numeric>
#include <queue>
#include <random>
#include <string>
#include <vector>
#include <atomic>
#ifdef _OPENMP
#include <omp.h>
#endif

using i32 = int32_t;
using i64 = int64_t;
using std::vector;
static const double INF = std::numeric_limits<double>::infinity();

static double now() {
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}

template <class T>
static vector<T> read_bin(const std::string& path) {
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    if (!f) { fprintf(stderr, "cannot open %s\n", path.c_str()); exit(1); }
    auto bytes = (size_t)f.tellg();
    f.seekg(0);
    vector<T> v(bytes / sizeof(T));
    f.read((char*)v.data(), bytes);
    return v;
}

// ------------------------------------------------------------------ graph

struct Graph {
    i64 n = 0, m = 0;
    vector<i64> off;          // directed CSR (original arcs)
    vector<i32> tgt;
    vector<double> x, y, z;   // z = elevation, if available
    vector<std::string> metrics;
    std::string dir;

    vector<i64> noff;         // undirected, de-duplicated
    vector<i32> nbr;
    vector<i64> rev;          // position of u in nbr list of v

    void load(const std::string& d) {
        dir = d;
        std::ifstream meta(d + "/meta.txt");
        std::string line;
        std::getline(meta, line); n = std::stoll(line);
        std::getline(meta, line); m = std::stoll(line);
        std::getline(meta, line);
        size_t p = 0;
        while (p < line.size()) {
            size_t q = line.find(' ', p);
            if (q == std::string::npos) q = line.size();
            if (q > p) metrics.push_back(line.substr(p, q - p));
            p = q + 1;
        }
        off = read_bin<i64>(d + "/off.i64");
        tgt = read_bin<i32>(d + "/tgt.i32");
        x = read_bin<double>(d + "/x.f64");
        y = read_bin<double>(d + "/y.f64");
        if (std::filesystem::exists(d + "/z.f64")) z = read_bin<double>(d + "/z.f64");
    }

    vector<double> weights(const std::string& metric) const {
        return read_bin<double>(dir + "/w_" + metric + ".f64");
    }

    void build_undirected() {
        vector<i64> deg(n + 1, 0);
        for (i64 u = 0; u < n; u++)
            for (i64 k = off[u]; k < off[u + 1]; k++) {
                i64 v = tgt[k];
                if (u != v) { deg[u + 1]++; deg[v + 1]++; }
            }
        for (i64 i = 0; i < n; i++) deg[i + 1] += deg[i];
        vector<i64> pos(deg.begin(), deg.end() - 1);
        vector<i32> tmp(deg[n]);
        for (i64 u = 0; u < n; u++)
            for (i64 k = off[u]; k < off[u + 1]; k++) {
                i64 v = tgt[k];
                if (u != v) { tmp[pos[u]++] = (i32)v; tmp[pos[v]++] = (i32)u; }
            }
        noff.assign(n + 1, 0);
        nbr.clear();
        nbr.reserve(deg[n]);
        for (i64 u = 0; u < n; u++) {
            std::sort(tmp.begin() + deg[u], tmp.begin() + deg[u + 1]);
            i32 last = -1;
            for (i64 k = deg[u]; k < deg[u + 1]; k++)
                if (tmp[k] != last) { nbr.push_back(tmp[k]); last = tmp[k]; }
            noff[u + 1] = (i64)nbr.size();
        }
        tmp.clear(); tmp.shrink_to_fit();
        rev.resize(nbr.size());
        for (i64 u = 0; u < n; u++)
            for (i64 k = noff[u]; k < noff[u + 1]; k++) {
                i64 v = nbr[k];
                rev[k] = std::lower_bound(nbr.begin() + noff[v], nbr.begin() + noff[v + 1],
                                          (i32)u) - nbr.begin();
            }
    }
};

// ------------------------------------------------- inertial flow bisection

struct Bisector {
    const Graph& g;
    vector<i64> member, level, ptr;
    vector<int8_t> role, intf, side;
    vector<i32> flow;
    i64 stamp = 0, lv = 0;

    explicit Bisector(const Graph& gg)
        : g(gg), member(gg.n, 0), level(2 * gg.n, 0), ptr(2 * gg.n, 0),
          role(gg.n, 0), intf(gg.n, 0), side(gg.n, 0), flow(gg.nbr.size(), 0) {}

    // min vertex cut between the first/last 25% along direction `dirn`
    bool split(vector<i64>& seg, int dirn, vector<i64>& sep, vector<i64>& A, vector<i64>& B) {
        const i64 m = (i64)seg.size();
        vector<i64> order(seg);
        auto key = [&](i64 v) {
            switch (dirn) {
                case 0: return g.x[v];
                case 1: return g.y[v];
                case 2: return g.x[v] + g.y[v];
                default: return g.x[v] - g.y[v];
            }
        };
        std::sort(order.begin(), order.end(), [&](i64 a, i64 b) { return key(a) < key(b); });
        const i64 k = std::max<i64>(1, m / 4);
        for (i64 i = 0; i < m; i++) role[order[i]] = 0;
        for (i64 i = 0; i < k; i++) { role[order[i]] = 1; role[order[m - 1 - i]] = 2; }
        for (i64 i = 0; i < k; i++) {
            i64 s = order[i];
            for (i64 j = g.noff[s]; j < g.noff[s + 1]; j++) {
                i64 v = g.nbr[j];
                if (member[v] == stamp && role[v] == 2) return false;  // adjacent source/sink
            }
        }
        for (i64 i = 0; i < m; i++) {
            i64 v = seg[i];
            intf[v] = 0;
            for (i64 j = g.noff[v]; j < g.noff[v + 1]; j++) flow[j] = 0;
        }
        vector<i64> q(2 * m), stack(2 * m + 2);
        auto lvl_of = [&](i64 node) { return level[node] & 0xFFFFFFFF; };
        auto seen = [&](i64 node) { return (level[node] >> 32) == lv; };

        while (true) {
            lv++;
            i64 qh = 0, qt = 0;
            for (i64 i = 0; i < k; i++) {
                i64 s = order[i];
                level[2 * s] = lv << 32; q[qt++] = 2 * s;
                level[2 * s + 1] = lv << 32; q[qt++] = 2 * s + 1;
            }
            bool reached = false;
            while (qh < qt) {
                i64 xn = q[qh++], v = xn >> 1, d = lvl_of(xn);
                if ((xn & 1) == 0 && role[v] == 2) { reached = true; continue; }
                if (xn & 1) {
                    if (role[v] == 0 && intf[v] == 1 && !seen(2 * v)) {
                        level[2 * v] = (lv << 32) | (d + 1); q[qt++] = 2 * v;
                    }
                    for (i64 j = g.noff[v]; j < g.noff[v + 1]; j++) {
                        i64 u = g.nbr[j];
                        if (member[u] == stamp && !seen(2 * u)) {
                            level[2 * u] = (lv << 32) | (d + 1); q[qt++] = 2 * u;
                        }
                    }
                } else {
                    if ((role[v] != 0 || intf[v] == 0) && !seen(2 * v + 1)) {
                        level[2 * v + 1] = (lv << 32) | (d + 1); q[qt++] = 2 * v + 1;
                    }
                    for (i64 j = g.noff[v]; j < g.noff[v + 1]; j++) {
                        i64 u = g.nbr[j];
                        if (member[u] == stamp && flow[j] < 0 && !seen(2 * u + 1)) {
                            level[2 * u + 1] = (lv << 32) | (d + 1); q[qt++] = 2 * u + 1;
                        }
                    }
                }
            }
            if (!reached) break;
            for (i64 i = 0; i < qt; i++) ptr[q[i]] = 0;
            for (i64 i = 0; i < k; i++) {
                const i64 start = 2 * order[i] + 1;
                while (true) {
                    i64 sp = 0;
                    stack[0] = start;
                    bool found = false;
                    while (sp >= 0) {
                        i64 xn = stack[sp], v = xn >> 1;
                        if ((xn & 1) == 0 && role[v] == 2) { found = true; break; }
                        const i64 d = lvl_of(xn), deg = g.noff[v + 1] - g.noff[v];
                        bool advanced = false;
                        while (ptr[xn] <= deg) {
                            const i64 a = ptr[xn];
                            bool ok;
                            i64 tn = 0;
                            if (a == 0) {
                                if (xn & 1) { ok = role[v] == 0 && intf[v] == 1; tn = 2 * v; }
                                else { ok = role[v] != 0 || intf[v] == 0; tn = 2 * v + 1; }
                            } else {
                                const i64 j = g.noff[v] + a - 1, u = g.nbr[j];
                                if (member[u] != stamp) ok = false;
                                else if (xn & 1) { ok = true; tn = 2 * u; }
                                else { ok = flow[j] < 0; tn = 2 * u + 1; }
                            }
                            if (ok && seen(tn) && lvl_of(tn) == d + 1) {
                                stack[++sp] = tn; advanced = true; break;
                            }
                            ptr[xn]++;
                        }
                        if (!advanced) {
                            level[xn] = 0;
                            if (--sp >= 0) ptr[stack[sp]]++;
                        }
                    }
                    if (!found) break;
                    for (i64 t = 0; t < sp; t++) {
                        const i64 a = stack[t], b = stack[t + 1], va = a >> 1, vb = b >> 1;
                        if (va == vb) {
                            if (role[va] == 0) intf[va] = (a & 1) ? 0 : 1;
                        } else {
                            const i64 j = g.noff[va] + ptr[a] - 1;
                            flow[j]++;
                            flow[g.rev[j]]--;
                        }
                    }
                }
            }
        }
        // residual reachability = min cut
        lv++;
        i64 qh = 0, qt = 0;
        for (i64 i = 0; i < k; i++) {
            i64 s = order[i];
            level[2 * s] = lv << 32; q[qt++] = 2 * s;
            level[2 * s + 1] = lv << 32; q[qt++] = 2 * s + 1;
        }
        while (qh < qt) {
            i64 xn = q[qh++], v = xn >> 1;
            if (xn & 1) {
                if (role[v] == 0 && intf[v] == 1 && !seen(2 * v)) { level[2 * v] = lv << 32; q[qt++] = 2 * v; }
                for (i64 j = g.noff[v]; j < g.noff[v + 1]; j++) {
                    i64 u = g.nbr[j];
                    if (member[u] == stamp && !seen(2 * u)) { level[2 * u] = lv << 32; q[qt++] = 2 * u; }
                }
            } else {
                if ((role[v] != 0 || intf[v] == 0) && !seen(2 * v + 1)) { level[2 * v + 1] = lv << 32; q[qt++] = 2 * v + 1; }
                for (i64 j = g.noff[v]; j < g.noff[v + 1]; j++) {
                    i64 u = g.nbr[j];
                    if (member[u] == stamp && flow[j] < 0 && !seen(2 * u + 1)) { level[2 * u + 1] = lv << 32; q[qt++] = 2 * u + 1; }
                }
            }
        }
        sep.clear(); A.clear(); B.clear();
        for (i64 i = 0; i < m; i++) {
            const i64 v = seg[i];
            const bool rin = seen(2 * v), rout = seen(2 * v + 1);
            if (rin && !rout) sep.push_back(v);
            else if (rout) A.push_back(v);
            else B.push_back(v);
        }
        return !A.empty() && !B.empty();
    }

    void geometric(vector<i64>& seg, vector<i64>& sep, vector<i64>& A, vector<i64>& B) {
        const i64 m = (i64)seg.size();
        double xmin = INF, xmax = -INF, ymin = INF, ymax = -INF;
        for (i64 v : seg) {
            xmin = std::min(xmin, g.x[v]); xmax = std::max(xmax, g.x[v]);
            ymin = std::min(ymin, g.y[v]); ymax = std::max(ymax, g.y[v]);
        }
        const bool byx = (xmax - xmin) >= (ymax - ymin);
        vector<i64> order(seg);
        std::sort(order.begin(), order.end(),
                  [&](i64 a, i64 b) { return byx ? g.x[a] < g.x[b] : g.y[a] < g.y[b]; });
        for (i64 i = 0; i < m; i++) side[order[i]] = (i < m / 2) ? 0 : 1;
        sep.clear(); A.clear(); B.clear();
        for (i64 i = 0; i < m; i++) {
            const i64 v = order[i];
            if (side[v] == 1) { B.push_back(v); continue; }
            bool boundary = false;
            for (i64 j = g.noff[v]; j < g.noff[v + 1] && !boundary; j++) {
                const i64 u = g.nbr[j];
                if (member[u] == stamp && side[u] == 1) boundary = true;
            }
            (boundary ? sep : A).push_back(v);
        }
    }
};

static vector<i32> nested_dissection(Graph& g, i64 leaf) {
    Bisector bi(g);
    vector<i64> buf(g.n);
    std::iota(buf.begin(), buf.end(), 0);
    vector<std::pair<i64, i64>> stack{{0, g.n}};
    vector<i64> seg, sep, A, B, sepb, Ab, Bb;
    vector<i64> deg(g.n);
    for (i64 v = 0; v < g.n; v++) deg[v] = g.noff[v + 1] - g.noff[v];
    while (!stack.empty()) {
        auto [lo, hi] = stack.back();
        stack.pop_back();
        const i64 m = hi - lo;
        if (m <= 0) continue;
        seg.assign(buf.begin() + lo, buf.begin() + hi);
        if (m <= leaf) {
            std::stable_sort(seg.begin(), seg.end(), [&](i64 a, i64 b) { return deg[a] < deg[b]; });
            std::copy(seg.begin(), seg.end(), buf.begin() + lo);
            continue;
        }
        bi.stamp++;
        for (i64 v : seg) bi.member[v] = bi.stamp;
        i64 best = -1, bestbal = 0;
        for (int d = 0; d < 4; d++) {
            if (!bi.split(seg, d, sep, A, B)) continue;
            const i64 bal = std::abs((i64)A.size() - (i64)B.size());
            if (best < 0 || (i64)sep.size() < best || ((i64)sep.size() == best && bal < bestbal)) {
                best = (i64)sep.size(); bestbal = bal; sepb = sep; Ab = A; Bb = B;
            }
        }
        if (best < 0) { bi.geometric(seg, sepb, Ab, Bb); }
        i64 p = lo;
        for (i64 v : Ab) buf[p++] = v;
        for (i64 v : Bb) buf[p++] = v;
        for (i64 v : sepb) buf[p++] = v;
        stack.push_back({lo, lo + (i64)Ab.size()});
        stack.push_back({lo + (i64)Ab.size(), lo + (i64)Ab.size() + (i64)Bb.size()});
    }
    vector<i32> rank(g.n);
    for (i64 i = 0; i < g.n; i++) rank[buf[i]] = (i32)i;
    return rank;
}

// ------------------------------------------------------------- contraction

struct CCH {
    i64 n = 0;
    vector<i32> rank;
    vector<i64> first;     // upward arcs per rank
    vector<i32> head;
    vector<i32> parent;
    vector<i64> e_arc;     // per original arc
    vector<int8_t> e_up;
    i64 depth = 0, triangles = 0;

    void symbolic(const Graph& g) {
        n = g.n;
        vector<i32> inv(n);
        for (i64 v = 0; v < n; v++) inv[rank[v]] = (i32)v;
        parent.assign(n, -1);
        vector<i32> child_head(n, -1), child_next(n, -1);
        vector<i64> marker(n, -1), tmp;
        first.assign(n + 1, 0);
        head.clear();
        for (i64 x = 0; x < n; x++) {
            const i64 v = inv[x];
            tmp.clear();
            for (i64 j = g.noff[v]; j < g.noff[v + 1]; j++) {
                const i64 r = rank[g.nbr[j]];
                if (r > x && marker[r] != x) { marker[r] = x; tmp.push_back(r); }
            }
            for (i32 c = child_head[x]; c >= 0; c = child_next[c])
                for (i64 a = first[c]; a < first[c + 1]; a++) {
                    const i64 r = head[a];
                    if (r != x && marker[r] != x) { marker[r] = x; tmp.push_back(r); }
                }
            std::sort(tmp.begin(), tmp.end());
            for (i64 r : tmp) head.push_back((i32)r);
            first[x + 1] = (i64)head.size();
            if (!tmp.empty()) {
                const i32 p = (i32)tmp[0];
                parent[x] = p;
                child_next[x] = child_head[p];
                child_head[p] = (i32)x;
            }
        }
        vector<i32> dep(n, 0);
        for (i64 v = n - 1; v >= 0; v--)
            if (parent[v] >= 0) { dep[v] = dep[parent[v]] + 1; depth = std::max<i64>(depth, dep[v]); }
        triangles = 0;
        for (i64 x = 0; x < n; x++) {
            const i64 d = first[x + 1] - first[x];
            triangles += d * (d - 1) / 2;
        }
    }

    i64 find_arc(i64 lo, i64 hi) const {
        return std::lower_bound(head.begin() + first[lo], head.begin() + first[lo + 1], (i32)hi)
               - head.begin();
    }

    void map_arcs(const Graph& g) {
        e_arc.assign(g.m, -1);
        e_up.assign(g.m, 0);
        for (i64 u = 0; u < g.n; u++) {
            const i64 ru = rank[u];
            for (i64 k = g.off[u]; k < g.off[u + 1]; k++) {
                const i64 rv = rank[g.tgt[k]];
                if (ru == rv) continue;
                e_arc[k] = ru < rv ? find_arc(ru, rv) : find_arc(rv, ru);
                e_up[k] = ru < rv ? 1 : 0;
            }
        }
    }

    void customize(const vector<double>& w, vector<double>& up, vector<double>& dn) const {
        up.assign(head.size(), INF);
        dn.assign(head.size(), INF);
        for (size_t k = 0; k < e_arc.size(); k++) {
            const i64 a = e_arc[k];
            if (a < 0) continue;
            double& slot = e_up[k] ? up[a] : dn[a];
            if (w[k] < slot) slot = w[k];
        }
        for (i64 x = 0; x < n; x++) {
            const i64 s = first[x], e = first[x + 1];
            for (i64 i = s; i < e; i++) {
                const i64 y = head[i];
                i64 j = first[y];
                const i64 je = first[y + 1];
                const double d_yx = dn[i], u_xy = up[i];
                for (i64 kk = i + 1; kk < e; kk++) {
                    const i64 z = head[kk];
                    while (j < je && head[j] < z) j++;
                    const double a1 = d_yx + up[kk];
                    if (a1 < up[j]) up[j] = a1;
                    const double a2 = dn[kk] + u_xy;
                    if (a2 < dn[j]) dn[j] = a2;
                }
            }
        }
    }

    // Level of a vertex: one more than the highest level among its lower
    // neighbours. Vertices of the same level have no dependencies between them,
    // and they only write arcs whose tail has a strictly higher level.
    vector<i32> levels(vector<vector<i32>>& buckets) const {
        vector<i32> lvl(n, 0);
        i32 maxlvl = 0;
        for (i64 x = 0; x < n; x++)
            for (i64 a = first[x]; a < first[x + 1]; a++) {
                const i64 y = head[a];
                if (lvl[x] + 1 > lvl[y]) lvl[y] = lvl[x] + 1;
                maxlvl = std::max(maxlvl, lvl[y]);
            }
        buckets.assign(maxlvl + 1, {});
        for (i64 x = 0; x < n; x++) buckets[lvl[x]].push_back((i32)x);
        return lvl;
    }

    static void atomic_min(double& slot, double v) {
        auto* p = reinterpret_cast<std::atomic<uint64_t>*>(&slot);
        uint64_t old = p->load(std::memory_order_relaxed);
        double cur;
        std::memcpy(&cur, &old, sizeof cur);
        while (v < cur) {
            uint64_t nv;
            std::memcpy(&nv, &v, sizeof nv);
            if (p->compare_exchange_weak(old, nv, std::memory_order_relaxed)) return;
            std::memcpy(&cur, &old, sizeof cur);
        }
    }

    // Same result as customize(), with the vertices of each level in parallel.
    void customize_parallel(const vector<double>& w, vector<double>& up, vector<double>& dn,
                            const vector<vector<i32>>& buckets) const {
        up.assign(head.size(), INF);
        dn.assign(head.size(), INF);
        for (size_t k = 0; k < e_arc.size(); k++) {
            const i64 a = e_arc[k];
            if (a < 0) continue;
            double& slot = e_up[k] ? up[a] : dn[a];
            if (w[k] < slot) slot = w[k];
        }
        for (const auto& bucket : buckets) {
            const int cnt = (int)bucket.size();
            // Small levels (the top of the hierarchy) are cheaper to run serially
            // than to synchronise; only wide levels are worth parallelising.
            const bool par = cnt >= 512;
#pragma omp parallel for schedule(dynamic, 64) if (par)
            for (int bi = 0; bi < cnt; bi++) {
                const i64 x = bucket[bi];
                const i64 s = first[x], e = first[x + 1];
                for (i64 i = s; i < e; i++) {
                    const i64 y = head[i];
                    i64 j = first[y];
                    const i64 je = first[y + 1];
                    const double d_yx = dn[i], u_xy = up[i];
                    for (i64 kk = i + 1; kk < e; kk++) {
                        const i64 z = head[kk];
                        while (j < je && head[j] < z) j++;
                        atomic_min(up[j], d_yx + up[kk]);
                        atomic_min(dn[j], dn[kk] + u_xy);
                    }
                }
            }
        }
    }

    i64 negative_cycle_arcs(const vector<double>& up, const vector<double>& dn) const {
        i64 c = 0;
        for (size_t a = 0; a < up.size(); a++)
            if (up[a] + dn[a] < -1e-6) c++;
        return c;
    }
};

// ------------------------------------------------------------------ queries

struct Query {
    const CCH& c;
    vector<double> df, dr;
    vector<i64> mark;
    i64 stamp = 0;
    explicit Query(const CCH& cc) : c(cc), df(cc.n, INF), dr(cc.n, INF), mark(cc.n, 0) {}

    // elimination-tree query (sweep in rank order; valid for negative weights)
    double run(i64 s, i64 t, const vector<double>& up, const vector<double>& dn) {
        const i64 rs = c.rank[s], rt = c.rank[t];
        stamp++;
        for (i64 x = rs; x >= 0; x = c.parent[x]) { df[x] = INF; mark[x] = stamp; }
        for (i64 x = rt; x >= 0; x = c.parent[x]) dr[x] = INF;
        df[rs] = 0.0;
        dr[rt] = 0.0;
        for (i64 x = rs; x >= 0; x = c.parent[x]) {
            const double dx = df[x];
            if (dx == INF) continue;
            for (i64 a = c.first[x]; a < c.first[x + 1]; a++) {
                const double v = dx + up[a];
                if (v < df[c.head[a]]) df[c.head[a]] = v;
            }
        }
        double best = INF;
        for (i64 x = rt; x >= 0; x = c.parent[x]) {
            const double dx = dr[x];
            if (dx == INF) continue;
            if (mark[x] == stamp && dx + df[x] < best) best = dx + df[x];
            for (i64 a = c.first[x]; a < c.first[x + 1]; a++) {
                const double v = dx + dn[a];
                if (v < dr[c.head[a]]) dr[c.head[a]] = v;
            }
        }
        return best;
    }
};

// Dijkstra on w + p[u] - p[v] (reference, and the query used by the shifted-CCH baseline)
static double dijkstra_potential(const Graph& g, const vector<double>& w, const vector<double>& p,
                                 i64 s, i64 t, vector<double>& dist, vector<char>& done, i64& scanned) {
    std::fill(dist.begin(), dist.end(), INF);
    std::fill(done.begin(), done.end(), 0);
    std::priority_queue<std::pair<double, i64>, vector<std::pair<double, i64>>, std::greater<>> pq;
    dist[s] = 0.0;
    pq.push({0.0, s});
    scanned = 0;
    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();
        if (done[u]) continue;
        done[u] = 1;
        scanned++;
        if (u == t) break;
        for (i64 k = g.off[u]; k < g.off[u + 1]; k++) {
            const i64 v = g.tgt[k];
            double rw = w[k] + p[u] - p[v];
            if (rw < 0) rw = 0;
            const double nd = d + rw;
            if (nd < dist[v]) { dist[v] = nd; pq.push({nd, v}); }
        }
    }
    return dist[t] == INF ? INF : dist[t] - p[s] + p[t];
}

// Johnson potential: FIFO (queue-based) Bellman-Ford from a virtual source with a
// 0-arc to every vertex. O(nm) worst case; `pops` counts queue removals.
// No negative-cycle detection: only call it on a conservative metric.
static vector<double> johnson_potential(const Graph& g, const vector<double>& w, i64& pops) {
    vector<double> p(g.n, 0.0);
    vector<char> inq(g.n, 1);
    std::deque<i64> q;
    for (i64 v = 0; v < g.n; v++) q.push_back(v);
    pops = 0;
    while (!q.empty()) {
        const i64 u = q.front();
        q.pop_front();
        pops++;
        inq[u] = 0;
        for (i64 k = g.off[u]; k < g.off[u + 1]; k++) {
            const i64 v = g.tgt[k];
            const double nd = p[u] + w[k];
            if (nd < p[v] - 1e-9) {
                p[v] = nd;
                if (!inq[v]) { inq[v] = 1; q.push_back(v); }
            }
        }
    }
    return p;
}

// The same virtual-source potential p(v) = min(0, min_u dist(u,v)), read off a
// customized CCH in two linear sweeps: every shortest u-v path has an up-down
// representative (Theorem 1), so an upward pass over all vertices in increasing
// rank followed by a downward pass in decreasing rank finds it. O(|E+|) time.
static vector<double> cch_potential(const CCH& c, const vector<double>& up, const vector<double>& dn) {
    vector<double> pr(c.n, 0.0);            // indexed by rank
    for (i64 x = 0; x < c.n; x++) {
        const double px = pr[x];
        for (i64 a = c.first[x]; a < c.first[x + 1]; a++) {
            const double v = px + up[a];
            if (v < pr[c.head[a]]) pr[c.head[a]] = v;
        }
    }
    for (i64 x = c.n - 1; x >= 0; x--) {
        double best = pr[x];
        for (i64 a = c.first[x]; a < c.first[x + 1]; a++) {
            const double v = pr[c.head[a]] + dn[a];
            if (v < best) best = v;
        }
        pr[x] = best;
    }
    vector<double> p(c.n);
    for (i64 v = 0; v < c.n; v++) p[v] = pr[c.rank[v]];
    return p;
}

// CCH query with a Dijkstra-based bidirectional search and stall-on-demand.
// Only valid on non-negative (e.g. potential-shifted) weights; used as the
// classical baseline for comparison.
struct ShiftedQuery {
    const CCH& c;
    vector<double> df, dr;
    vector<i64> stf, str_;
    i64 stamp = 0;
    explicit ShiftedQuery(const CCH& cc) : c(cc), df(cc.n, INF), dr(cc.n, INF),
                                           stf(cc.n, 0), str_(cc.n, 0) {}
    double run(i64 s, i64 t, const vector<double>& up, const vector<double>& dn, bool stall,
               i64& settled) {
        stamp++;
        std::priority_queue<std::pair<double, i64>, vector<std::pair<double, i64>>, std::greater<>> qf, qr;
        const i64 rs = c.rank[s], rt = c.rank[t];
        df[rs] = 0.0; stf[rs] = stamp; qf.push({0.0, rs});
        dr[rt] = 0.0; str_[rt] = stamp; qr.push({0.0, rt});
        double mu = INF;
        settled = 0;
        while (!qf.empty() || !qr.empty()) {
            const double tf = qf.empty() ? INF : qf.top().first;
            const double tr = qr.empty() ? INF : qr.top().first;
            if (std::min(tf, tr) >= mu) break;
            const bool fwd = tf <= tr;
            auto& q = fwd ? qf : qr;
            auto& d = fwd ? df : dr;
            auto& st = fwd ? stf : str_;
            auto& od = fwd ? dr : df;
            auto& ost = fwd ? str_ : stf;
            const auto& wf = fwd ? up : dn;      // forward uses upward weights
            const auto& wb = fwd ? dn : up;      // for stalling
            auto [dx, x] = q.top();
            q.pop();
            if (st[x] != stamp || dx > d[x]) continue;
            if (stall) {                          // stall-on-demand
                bool stalled = false;
                for (i64 a = c.first[x]; a < c.first[x + 1] && !stalled; a++) {
                    const i64 y = c.head[a];
                    if (st[y] == stamp && d[y] + wb[a] < dx - 1e-12) stalled = true;
                }
                if (stalled) continue;
            }
            settled++;
            if (ost[x] == stamp && dx + od[x] < mu) mu = dx + od[x];
            for (i64 a = c.first[x]; a < c.first[x + 1]; a++) {
                const i64 y = c.head[a];
                const double nd = dx + wf[a];
                if (st[y] != stamp || nd < d[y]) { d[y] = nd; st[y] = stamp; q.push({nd, y}); }
            }
        }
        return mu;
    }
};

// ------------------------------------------------------------------- main

// Repeated measurements: median and range, in milliseconds.
struct Timing {
    vector<double> ms;
    void add(double seconds) { ms.push_back(1000.0 * seconds); }
    double med() const {
        vector<double> v = ms;
        std::sort(v.begin(), v.end());
        const size_t k = v.size();
        return k % 2 ? v[k / 2] : 0.5 * (v[k / 2 - 1] + v[k / 2]);
    }
    std::string str() const {
        char b[160];
        snprintf(b, sizeof b, "%.4f ms (median of %zu, range %.4f-%.4f)", med(), ms.size(),
                 *std::min_element(ms.begin(), ms.end()), *std::max_element(ms.begin(), ms.end()));
        return b;
    }
};

static bool same(double a, double b) {
    return a == b || std::abs(a - b) <= 1e-6 * std::max(1.0, std::abs(b));
}

int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: cch <graph-dir> [--order f] [--queries n] [--check n] [--metric m] "
                        "[--shifted-cch] [--parallel threads] [--repeat r]\n");
        return 1;
    }
    std::string dir = argv[1], order_file, only_metric;
    i64 nq = 1000, ncheck = 20;
    bool shifted_cch = false;
    int threads = 0, reps = 5;
    for (int i = 2; i < argc; i++) {
        std::string a = argv[i];
        if (a == "--order" && i + 1 < argc) order_file = argv[++i];
        else if (a == "--queries" && i + 1 < argc) nq = std::stoll(argv[++i]);
        else if (a == "--check" && i + 1 < argc) ncheck = std::stoll(argv[++i]);
        else if (a == "--metric" && i + 1 < argc) only_metric = argv[++i];
        else if (a == "--shifted-cch") shifted_cch = true;
        else if (a == "--parallel" && i + 1 < argc) threads = std::stoi(argv[++i]);
        else if (a == "--repeat" && i + 1 < argc) reps = std::max(1, std::stoi(argv[++i]));
    }

    Graph g;
    double t0 = now();
    g.load(dir);
    printf("graph: n=%lld m=%lld metrics=%zu (%.1fs)\n", (long long)g.n, (long long)g.m,
           g.metrics.size(), now() - t0);
    t0 = now();
    g.build_undirected();
    printf("undirected: %zu slots (%.1fs)\n", g.nbr.size(), now() - t0);

    CCH c;
    if (!order_file.empty() && std::filesystem::exists(order_file)) {
        c.rank = read_bin<i32>(order_file);
        printf("order: loaded from %s\n", order_file.c_str());
    } else {
        t0 = now();
        c.rank = nested_dissection(g, 64);
        printf("order: %.1fs\n", now() - t0);
        if (!order_file.empty()) {
            std::ofstream f(order_file, std::ios::binary);
            f.write((const char*)c.rank.data(), c.rank.size() * sizeof(i32));
        }
    }
    t0 = now();
    c.symbolic(g);
    const double t_sym = now() - t0;
    t0 = now();
    c.map_arcs(g);
    printf("contraction: arcs=%zu triangles=%lld depth=%lld (symbolic %.1fs, mapping %.1fs)\n",
           c.head.size(), (long long)c.triangles, (long long)c.depth, t_sym, now() - t0);
    printf("repetitions per timing: %d; queries: %lld random pairs; verified: %lld\n", reps,
           (long long)nq, (long long)ncheck);

    std::mt19937_64 rng(7);
    vector<std::pair<i64, i64>> pairs;
    for (i64 i = 0; i < nq; i++)
        pairs.push_back({(i64)(rng() % g.n), (i64)(rng() % g.n)});

    vector<double> up, dn, dist(g.n);
    vector<char> done(g.n);
    Query q(c);
    ShiftedQuery sq(c);
    for (const auto& metric : g.metrics) {
        if (!only_metric.empty() && metric != only_metric) continue;
        const vector<double> w = g.weights(metric);
        const char* M = metric.c_str();

        i64 neg_arcs = 0, neg_loops = 0;
        for (i64 u = 0; u < g.n; u++)
            for (i64 k = g.off[u]; k < g.off[u + 1]; k++) {
                if (w[k] < 0) neg_arcs++;
                if (g.tgt[k] == u && w[k] < 0) neg_loops++;   // loops never enter G+
            }
        printf("[%s] negative arcs: %lld (%.1f%%), negative self-loops: %lld\n", M,
               (long long)neg_arcs, 100.0 * neg_arcs / g.m, (long long)neg_loops);

        Timing t_cust;
        for (int r = 0; r < reps; r++) {
            t0 = now();
            c.customize(w, up, dn);
            t_cust.add(now() - t0);
        }
        Timing t_neg;
        i64 neg = 0;
        for (int r = 0; r < reps; r++) {
            t0 = now();
            neg = c.negative_cycle_arcs(up, dn);
            t_neg.add(now() - t0);
        }
        printf("[%s] customization (no potential): %s\n", M, t_cust.str().c_str());
        printf("[%s] negative-cycle scan: %s, %lld arcs with a negative 2-cycle\n", M,
               t_neg.str().c_str(), (long long)neg);

        // Level-parallel customization, run last so that its multi-threaded load
        // does not affect the single-threaded timings above.
        auto parallel_check = [&]() {
            if (threads <= 0) return;
#ifdef _OPENMP
            omp_set_num_threads(threads);
#endif
            vector<vector<i32>> buckets;
            c.levels(buckets);
            vector<double> pup, pdn;
            Timing t_par;
            for (int r = 0; r < reps; r++) {
                t0 = now();
                c.customize_parallel(w, pup, pdn, buckets);
                t_par.add(now() - t0);
            }
            i64 diff = 0;   // exact comparison (INF == INF holds; no NaN can occur)
            for (size_t a = 0; a < up.size(); a++)
                if (!(up[a] == pup[a] && dn[a] == pdn[a])) diff++;
            printf("[%s] parallel customization (%d threads, %zu levels): %s, speed-up %.2fx, "
                   "arcs differing from serial (exact comparison): %lld\n",
                   M, threads, buckets.size(), t_par.str().c_str(), t_cust.med() / t_par.med(),
                   (long long)diff);
        };

        Timing t_q;
        i64 reach = 0;
        for (auto [s, t] : pairs) q.run(s, t, up, dn);   // warm-up, not timed
        for (int r = 0; r < reps; r++) {
            reach = 0;
            t0 = now();
            for (auto [s, t] : pairs) reach += q.run(s, t, up, dn) == INF ? 0 : 1;
            t_q.add((now() - t0) / (double)pairs.size());
        }
        printf("[%s] elimination-tree query (no potential): %s per query, %lld/%zu reachable\n", M,
               t_q.str().c_str(), (long long)reach, pairs.size());

        Timing t_pc;
        vector<double> pc;
        for (int r = 0; r < reps; r++) {
            t0 = now();
            pc = cch_potential(c, up, dn);
            t_pc.add(now() - t0);
        }
        double worst_pc = 0;
        for (i64 u = 0; u < g.n; u++)
            for (i64 k = g.off[u]; k < g.off[u + 1]; k++)
                if (g.tgt[k] != u) worst_pc = std::min(worst_pc, w[k] + pc[u] - pc[g.tgt[k]]);
        printf("[%s] potential read off the CCH (two sweeps): %s, min reduced weight %.3g\n", M,
               t_pc.str().c_str(), worst_pc);

        if (ncheck <= 0) { parallel_check(); continue; }

        Timing t_pot;
        vector<double> p;
        i64 pops = 0;
        for (int r = 0; r < reps; r++) {
            t0 = now();
            p = johnson_potential(g, w, pops);
            t_pot.add(now() - t0);
        }
        double maxdiff = 0;
        for (i64 v = 0; v < g.n; v++) maxdiff = std::max(maxdiff, std::abs(p[v] - pc[v]));
        printf("[%s] Johnson potential (FIFO Bellman-Ford): %s, %.2f queue pops per vertex, "
               "max |p_BF - p_CCH| = %.3g\n", M, t_pot.str().c_str(), (double)pops / g.n, maxdiff);

        i64 ok = 0, scanned = 0;
        double t_ref = 0;
        for (i64 i = 0; i < ncheck && i < (i64)pairs.size(); i++) {
            const auto [s, t] = pairs[i];
            const double got = q.run(s, t, up, dn);
            const double t1 = now();
            const double want = dijkstra_potential(g, w, p, s, t, dist, done, scanned);
            t_ref += now() - t1;
            if (same(got, want)) ok++;
        }
        const i64 nchk = std::min<i64>(ncheck, (i64)pairs.size());
        printf("[%s] verified against Dijkstra with Johnson potential: %lld/%lld correct "
               "(reference %.1f ms per query)\n", M, (long long)ok, (long long)nchk,
               1000 * t_ref / (double)nchk);

        if (!shifted_cch) { parallel_check(); continue; }

        // Classical pipelines on the same hierarchy: shift by a potential, customize the
        // shifted (non-negative) metric, and query it with the elimination-tree query or
        // with the Dijkstra-based bidirectional CCH query (with and without stall-on-demand).
        auto baseline = [&](const char* name, const vector<double>& pot, bool clamp) {
            vector<double> ws(w.size());
            for (i64 u = 0; u < g.n; u++)
                for (i64 k = g.off[u]; k < g.off[u + 1]; k++) {
                    const double r = w[k] + pot[u] - pot[g.tgt[k]];
                    ws[k] = clamp ? std::max(0.0, r) : r;
                }
            vector<double> sup, sdn;
            Timing t_cs;
            for (int r = 0; r < reps; r++) {
                t0 = now();
                c.customize(ws, sup, sdn);
                t_cs.add(now() - t0);
            }
            printf("[%s] %s-shifted: customization %s\n", M, name, t_cs.str().c_str());
            i64 agree = 0;
            for (i64 i = 0; i < nchk; i++) {
                const auto [s, t] = pairs[i];
                const double got = q.run(s, t, sup, sdn);
                const double back = got == INF ? INF : got - pot[s] + pot[t];
                if (same(back, q.run(s, t, up, dn))) agree++;
            }
            Timing t_et;
            for (auto [s, t] : pairs) q.run(s, t, sup, sdn);   // warm-up, not timed
            for (int r = 0; r < reps; r++) {
                t0 = now();
                for (auto [s, t] : pairs) q.run(s, t, sup, sdn);
                t_et.add((now() - t0) / (double)pairs.size());
            }
            printf("[%s] %s-shifted + elimination-tree query: %s per query, agrees %lld/%lld\n", M,
                   name, t_et.str().c_str(), (long long)agree, (long long)nchk);
            for (int stall = 0; stall < 2; stall++) {
                i64 settled_total = 0, okk = 0;
                Timing t_dq;
                for (auto [s, t] : pairs) {                   // warm-up, not timed
                    i64 settled = 0;
                    sq.run(s, t, sup, sdn, stall == 1, settled);
                }
                for (int r = 0; r < reps; r++) {
                    settled_total = 0;
                    t0 = now();
                    for (auto [s, t] : pairs) {
                        i64 settled = 0;
                        sq.run(s, t, sup, sdn, stall == 1, settled);
                        settled_total += settled;
                    }
                    t_dq.add((now() - t0) / (double)pairs.size());
                }
                for (i64 i = 0; i < nchk; i++) {
                    const auto [s, t] = pairs[i];
                    i64 settled = 0;
                    const double got = sq.run(s, t, sup, sdn, stall == 1, settled);
                    const double back = got == INF ? INF : got - pot[s] + pot[t];
                    if (same(back, q.run(s, t, up, dn))) okk++;
                }
                printf("[%s] %s-shifted + Dijkstra-based CCH query%s: %s per query, "
                       "%lld settled/query, agrees %lld/%lld\n", M, name,
                       stall ? " + stall-on-demand" : "", t_dq.str().c_str(),
                       (long long)(settled_total / (i64)pairs.size()), (long long)okk, (long long)nchk);
            }
        };

        // height-induced potential (free; the classical EV approach) where elevations exist
        vector<double> zm = g.z;
        if (std::filesystem::exists(dir + "/z_" + metric + ".f64"))
            zm = read_bin<double>(dir + "/z_" + metric + ".f64");
        else if (metric != "ev") zm.clear();
        if (!zm.empty()) {
            const double BETA_DOWN = 0.6;
            vector<double> ph(g.n);
            for (i64 v = 0; v < g.n; v++) ph[v] = BETA_DOWN * zm[v];
            double worst = 0;
            for (i64 u = 0; u < g.n; u++)
                for (i64 k = g.off[u]; k < g.off[u + 1]; k++)
                    worst = std::min(worst, w[k] + ph[u] - ph[g.tgt[k]]);
            printf("[%s] height-induced potential: free, min reduced weight %.3g\n", M, worst);
            baseline("height", ph, true);
        }
        baseline("Johnson", p, false);
        parallel_check();
    }
    return 0;
}
