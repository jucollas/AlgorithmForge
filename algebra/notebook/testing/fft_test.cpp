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
const int mod=998244353;
#include"../../modulo_int.cpp"
#include"../fft_icpc.cpp"
vector<mint> gen_vec(int n){
    vector<mint>vc(n);for(mint&ac:vc)ac=rng_64();
    return vc;
}vector<mint> brute(const vector<mint>&A,const vector<mint>&B){
    vector<mint>C(A.size()+B.size()-1,0);
    rep(i,0,A.size())rep(j,0,B.size())C[i+j]+=A[i]*B[j];
    return C;
}
void solve() {
    int n=1<<8,am=100,m=n+33;
    for(int tm=0;tm<am;++tm){
        vector<mint> A=gen_vec(n),B=gen_vec(n);
        assert(convolution(A,B)==brute(A,B));
        vector<mint> Ai=inverse(A,m);
        vector<mint> one=convolution(A,Ai);one.resize(m);
        assert(one[0]==1);rep(i,1,m)assert(one[i]==0);
    } cout << "FFT TEST: todo e bem" << endl;
} int32_t main(){
	ios_base::sync_with_stdio(false);
    cin.tie(NULL);
	cout << setprecision(12) << fixed;
#ifdef OSVARP
    auto start = chrono::high_resolution_clock::now();
#endif

    int t = 1;
    // cin >> t; 
    ++t; while ( --t ) solve();

#ifdef OSVARP
    auto end = chrono::high_resolution_clock::now();
    struct rusage usage; getrusage(RUSAGE_SELF, &usage);
    cerr << "\n<Execution time: "
        << chrono::duration_cast<chrono::milliseconds>(end - start).count()
        << " ms>\n<Memory used: "
        << usage.ru_maxrss << " kilobytes>" << endl;
#endif
    return 0; }