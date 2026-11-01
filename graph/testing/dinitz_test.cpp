// Author: Claude
// Tests for graph/dinitz.cpp
// Build: g++ -std=c++17 -O2 -fsanitize=address,undefined -o /tmp/dt dinitz_test.cpp && ./dt
// Checks, on many random graphs:
//   * max_flow matches an independent Edmonds-Karp
//   * the flow left in 'edges' is feasible (capacity bounds + conservation)
//   * the residual capacities c1/c2 are consistent with f
//   * the source side of the residual graph is a cut with capacity == flow
#include <bits/stdc++.h>
#include <cassert>
using namespace std;
#define all(vec) vec.begin(), vec.end()
#include "../dinitz.cpp"

mt19937_64 rng(12345);
long long rnd(long long lo, long long hi) { return lo + rng() % (hi - lo + 1); }

struct RawEdge { int u, v; long long c1, c2; };

// Edmonds-Karp on an adjacency matrix (parallel edges are merged)
long long reference(int n, const vector<RawEdge>& es, int s, int t) {
	vector<vector<long long>> cap(n, vector<long long>(n, 0));
	for (auto& e : es) if (e.u != e.v) cap[e.u][e.v] += e.c1, cap[e.v][e.u] += e.c2;
	long long total = 0;
	while (true) {
		vector<int> par(n, -1); par[s] = s; queue<int> q; q.push(s);
		while (!q.empty() && par[t] == -1) {
			int u = q.front(); q.pop();
			for (int v = 0; v < n; ++v) if (par[v] == -1 && cap[u][v] > 0) par[v] = u, q.push(v);
		}
		if (par[t] == -1) return total;
		long long b = LLONG_MAX;
		for (int v = t; v != s; v = par[v]) b = min(b, cap[par[v]][v]);
		for (int v = t; v != s; v = par[v]) cap[par[v]][v] -= b, cap[v][par[v]] += b;
		total += b;
	}
}

template <typename ftype, ftype INF>
void check(int n, const vector<RawEdge>& es, int s, int t, const char* tag) {
	Dinitz<ftype, INF> d(n); d.source = s; d.sink = t;
	for (auto& e : es) d.add_edge(e.u, e.v, ftype(e.c1), ftype(e.c2));
	ftype got = d.max_flow();
	long long want = reference(n, es, s, t);
	auto fail = [&](const string& why) {
		printf("FAIL [%s] %s\n n=%d s=%d t=%d got=%lld want=%lld\n", tag, why.c_str(), n, s, t,
		       (long long)got, want);
		for (auto& e : es) printf("  %d %d %lld %lld\n", e.u, e.v, e.c1, e.c2);
		exit(1);
	};
	if ((long long)got != want) fail("wrong max flow value");
	vector<long long> bal(n, 0);
	for (size_t i = 0; i < es.size(); ++i) {
		auto& e = d.edges[i]; auto& r = es[i];
		if (e.u != e.v) {
			if (e.f > r.c1 || -e.f > r.c2) fail("capacity violated on edge " + to_string(i));
			if (e.c1 != r.c1 - e.f || e.c2 != r.c2 + e.f) fail("residual inconsistent on edge " + to_string(i));
		}
		bal[e.u] -= e.f; bal[e.v] += e.f;
	}
	for (int v = 0; v < n; ++v) if (v != s && v != t && bal[v] != 0) fail("conservation violated at " + to_string(v));
	if (bal[t] != want || bal[s] != -want) fail("net flow at source/sink != reported value");
	// min cut from the residual graph
	vector<bool> vis(n, false); vector<int> st = {s}; vis[s] = true;
	while (!st.empty()) {
		int u = st.back(); st.pop_back();
		for (int id : d.graph[u]) { int v = d.edges[id].to(u); if (!vis[v] && d.edges[id].cap(u) > 0) vis[v] = true, st.push_back(v); }
	}
	if (vis[t]) fail("sink still reachable in residual graph");
	long long cut = 0;
	for (auto& e : es) { if (vis[e.u] && !vis[e.v]) cut += e.c1; if (vis[e.v] && !vis[e.u]) cut += e.c2; }
	if (cut != want) fail("residual cut capacity != flow");
}

vector<RawEdge> random_graph(int n, int m, long long maxc, bool bidir, bool zeros) {
	vector<RawEdge> es;
	for (int i = 0; i < m; ++i) {
		RawEdge e{int(rnd(0, n - 1)), int(rnd(0, n - 1)), rnd(zeros ? 0 : 1, maxc), bidir ? rnd(0, maxc) : 0};
		es.push_back(e);
	}
	return es;
}

int main() {
	constexpr int IINF = 1e9 + 33;
	constexpr long long LINF = (long long)4e18;
	int cnt = 0;
	// fixed hand-made cases
	check<int, IINF>(2, {}, 0, 1, "no edges");
	check<int, IINF>(2, {{0, 1, 5, 0}}, 0, 1, "single edge");
	check<int, IINF>(2, {{1, 0, 5, 0}}, 0, 1, "single edge wrong way");
	check<int, IINF>(2, {{1, 0, 0, 7}}, 0, 1, "single edge via c2");
	check<int, IINF>(3, {{0, 0, 9, 9}, {0, 1, 3, 0}, {1, 1, 4, 0}, {1, 2, 2, 0}}, 0, 2, "self loops");
	check<int, IINF>(2, {{0, 1, 3, 0}, {0, 1, 4, 0}, {1, 0, 0, 5}}, 0, 1, "parallel edges");
	check<int, IINF>(4, {{0, 1, 1, 0}, {0, 2, 1, 0}, {1, 3, 1, 0}, {2, 3, 1, 0}, {1, 2, 1, 1}}, 0, 3, "diamond");
	check<int, IINF>(4, {{0, 1, 0, 0}, {1, 2, 0, 0}, {2, 3, 0, 0}}, 0, 3, "all zero");
	check<int, IINF>(5, {{0, 1, 5, 0}, {1, 2, 5, 0}}, 0, 4, "isolated sink");
	check<long long, LINF>(3, {{0, 1, (long long)1e18, 0}, {1, 2, (long long)1e18, 0}, {0, 2, (long long)1e18, 0}}, 0, 2, "huge capacities");
	check<int, IINF>(3, {{0, 1, (1 << 30), 0}, {1, 2, (1 << 30) - 1, 0}}, 0, 2, "near int limit");
	// random cases
	for (int it = 0; it < 30000; ++it) {
		int n = rnd(2, 9), m = rnd(0, 25);
		int s = rnd(0, n - 1), t = rnd(0, n - 2); if (t >= s) ++t;
		long long maxc = vector<long long>{1, 2, 3, 10, 1000, 1000000}[rnd(0, 5)];
		auto es = random_graph(n, m, maxc, rnd(0, 1), rnd(0, 1));
		check<int, IINF>(n, es, s, t, "random small int"); ++cnt;
	}
	for (int it = 0; it < 2000; ++it) {
		int n = rnd(2, 60), m = rnd(0, 400);
		int s = rnd(0, n - 1), t = rnd(0, n - 2); if (t >= s) ++t;
		auto es = random_graph(n, m, (long long)rnd(1, 1e15), rnd(0, 1), rnd(0, 1));
		check<long long, LINF>(n, es, s, t, "random medium long long"); ++cnt;
	}
	// a layered graph (lots of equal-length paths) and a long path (deep recursion)
	{
		int L = 30, W = 30, n = L * W + 2; vector<RawEdge> es;
		for (int w = 0; w < W; ++w) es.push_back({0, 2 + w, rnd(1, 100), 0}), es.push_back({2 + (L - 1) * W + w, 1, rnd(1, 100), 0});
		for (int l = 0; l + 1 < L; ++l) for (int a = 0; a < W; ++a) for (int k = 0; k < 3; ++k)
			es.push_back({2 + l * W + a, 2 + (l + 1) * W + int(rnd(0, W - 1)), rnd(1, 100), 0});
		check<int, IINF>(n, es, 0, 1, "layered");
		int P = 20000; vector<RawEdge> path;
		for (int i = 0; i + 1 < P; ++i) path.push_back({i, i + 1, rnd(1, 1e9), rnd(0, 1e9)});
		check<long long, LINF>(P, path, 0, P - 1, "long path");
	}
	printf("OK (%d random cases + fixed cases)\n", cnt);
}
