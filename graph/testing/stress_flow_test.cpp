/*
 ________
|    ___ |
|  ,',.(`|
| :  `'  |
| :) _  (|
|  `:_)_,|
|________|

Autor: Oscar Vargas Pabon
Fecha: 

*/
//#pragma GCC optimize("O3")
//#define NDEBUG
#include <bits/stdc++.h>
#include <cassert>

typedef long long lint;
// __uint128_t, __int128_t, uint64_t, int64_t, uint32_t,int32_t
using namespace std;
#ifdef OSVARP
    #include<sys/resource.h>
#else
    #define cerr for(;false;) cerr
#endif
template <typename t1,typename t2> std::ostream &operator<<(ostream &os, const std::pair<t1,t2>&pr){return os<<"("<<pr.first<<";"<<pr.second<<")";};
template <typename T> std::ostream &operator<<(ostream &os, const std::vector<T>&vc){for(auto ac:vc)os<<ac<<' ';return os;};
template <typename T> std::ostream &operator<<(ostream &os, const std::set<   T>&vc){for(auto ac:vc)os<<ac<<' ';return os;};
template <typename t1,typename t2> std::ostream &operator<<(ostream &os, const std::map<t1,t2>&vc){for(auto ac:vc)os<<ac<<' ';return os;};
#define debug(args...) { string _s = #args; replace(_s.begin(), _s.end(), ',', ' '); stringstream _ss(_s); istream_iterator<string> _it(_ss); raw_debug(_it, args);}
void raw_debug(istream_iterator<string> it) {cerr<<endl;assert(it==it);}
template<typename T, typename... Args>
void raw_debug(istream_iterator<string> it, T a, Args... args) { cerr <<"<"<< *it << "->" << a << "> "; raw_debug(++it, args...); }
#define adebug(ar,n) {cerr<<'['<<#ar<<']';for(int my_imp_ind=0;my_imp_ind<n;++my_imp_ind)cerr << ' ' << ar[my_imp_ind]; cerr << endl;}

#define rep(i,strt,end) for(int i = strt ; i !=int(end) ; (int(strt)<int(end))?++i:--i )
#define rall(vec) vec.rbegin(), vec.rend()
#define all(vec) vec.begin(), vec.end()
#define eb emplace_back
#define pb push_back
#define pob pop_back
#define pf push_front
#define pof pop_front

std::mt19937_64 rng_64( std::chrono::steady_clock::now().time_since_epoch().count() );
constexpr int ilog2( int num ) { return 8*sizeof(int) - __builtin_clz( num ) - 1; }
template<typename tpow,typename texp=int64_t> constexpr tpow mpow(tpow x,uint64_t e,tpow m){tpow res=1;while(e){if(e&1)res=(texp(res)*x)%m;e>>=1;x=(texp(x)*x)%m;}return res;}

inline int gen(int limit){return (limit+rng_64()%limit)%limit;}
inline bool prob(int x,int y){return gen(y)<x;}//(x/y)
vector<pair<int,int>> erdos_renyi(int n,int x,int y){
    vector<pair<int,int>> edg;
    rep(i,0,n)rep(j,0,n)if(prob(x,y))edg.emplace_back(i,j);
    return edg;
}

#include"../dinitz.cpp"
namespace ek{
#include"../old_impl/edmondsKarps.cpp"
};
const int inf=1e9+33, tope=1e6;
void solve(int tm) {if(tm%10==0 )debug(tm);
    const int n=100,x=67,y=100;
    vector<pair<int,int>> redg=erdos_renyi(n,x,y);
    Dinitz<int,inf> fl1(n);
    int src=gen(n),snk=gen(n);while(snk==src)snk=gen(n);
    ek::graph.clear();ek::graph.resize(n);ek::edges.clear();
    ek::source=fl1.source=src;ek::sink=fl1.sink=snk;
    for(auto[u,v]:redg){
        const int c=gen(tope);
        fl1.add_edge(u,v,c);
        ek::addEdge(u,v,c);
    }int r1=fl1.max_flow(),r2=ek::maxFlow();
    assert(r1==r2);
	
} int32_t main(){
	ios_base::sync_with_stdio(false);
    cin.tie(NULL);
	cout << setprecision(12) << fixed;
#ifdef OSVARP
    auto start = chrono::high_resolution_clock::now();
#endif

    int t = 1e3;
    // cin >> t; 
    ++t; while ( --t ) solve(t);

#ifdef OSVARP
    auto end = chrono::high_resolution_clock::now();
    struct rusage usage; getrusage(RUSAGE_SELF, &usage);
    cerr << "\n<Execution time: "
        << chrono::duration_cast<chrono::milliseconds>(end - start).count()
        << " ms>\n<Memory used: "
        << usage.ru_maxrss << " kilobytes>" << endl;
#endif
    return 0; }