// Author: Claude
// Timing benchmark for graph/dinitz.cpp
// Build: g++ -std=c++17 -O2 -o /tmp/db dinitz_bench.cpp && /tmp/db [--slow]
// --slow also runs the Dinitz killer from https://codeforces.com/blog/entry/52714?#comment-397861
// (~60s with the current implementation: 1998 BFS phases on a 4M-edge graph)
#include <bits/stdc++.h>
#include <cassert>
using namespace std;
#define all(vec) vec.begin(), vec.end()
#include "../dinitz.cpp"
mt19937 rng(7);
struct E { int u, v; long long c; };
void run(const char* name, int n, const vector<E>& es, int s, int t) {
	Dinitz<long long, (long long)4e18> d(n); d.source = s; d.sink = t;
	for (auto& e : es) d.add_edge(e.u, e.v, e.c);
	auto st = chrono::steady_clock::now(); long long f = d.max_flow();
	printf("%-22s flow=%-16lld %6lld ms\n", name, f,
	       (long long)chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - st).count());
}
int main(int argc, char** argv) {
	{ int n = 2000; vector<E> v; for (int i = 0; i < 200000; i++) v.push_back({int(rng() % n), int(rng() % n), (long long)(rng() % 1000000000) + 1});
	  run("random sparse 1e9", n, v, 0, 1); }
	{ int n = 500; vector<E> v; for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (i != j && rng() % 2) v.push_back({i, j, (long long)(rng() % 1000000) + 1});
	  run("random dense 1e6", n, v, 0, n - 1); }
	{ int n = 5000; vector<E> v; for (int i = 0; i < n; i++) for (int k = 0; k < 20; k++)
	    v.push_back({i, int(rng() % n), 1 + (rng() % 3 == 0 ? (long long)(rng() % 1000000000) : (long long)(rng() % 3))});
	  run("mixed tiny/huge caps", n, v, 0, n - 1); }
	{ int L = 100, W = 100, n = L * W + 2; vector<E> v;
	  for (int w = 0; w < W; w++) v.push_back({0, 2 + w, (long long)4e15}), v.push_back({2 + (L - 1) * W + w, 1, (long long)4e15});
	  for (int l = 0; l + 1 < L; l++) for (int a = 0; a < W; a++) for (int k = 0; k < 5; k++)
	    v.push_back({2 + l * W + a, 2 + (l + 1) * W + int(rng() % W), (long long)(rng() % 1000000000000LL) + 1});
	  for (int k = 0; k < 50000; k++) v.push_back({int(rng() % (n - 2)) + 2, int(rng() % (n - 2)) + 2, (long long)(rng() % 1000) + 1});
	  run("layered + shortcuts", n, v, 0, 1); }
	if (argc > 1 && string(argv[1]) == "--slow") {
		vector<E> v;
		for (int i = 1; i <= 1998; i++) {
			v.push_back({i, i + 1, 100000000}); v.push_back({i, 2000, 10000});
			for (int j = i + 3; j <= 1998; j++) v.push_back({i, j, 1}), v.push_back({j, i, 1});
		}
		run("cf-52714 killer", 2001, v, 1, 2000);
	} // https://codeforces.com/blog/entry/52714
}
