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
//#pragma GCC optimize("Ofast")
//#define NDEBUG
#include <bits/stdc++.h>
#include <cassert>

typedef long long lint;
typedef __int128_t int128;
using namespace std;
#ifdef OSVARP
    #include<sys/resource.h>
#else
    #define cerr for(;false;) cerr
#endif
#define debug(args...) { string _s = #args; replace(_s.begin(), _s.end(), ',', ' '); stringstream _ss(_s); istream_iterator<string> _it(_ss); raw_debug(_it, args);}
void raw_debug(istream_iterator<string> it) {cerr<<endl;assert(it==it);}
template<typename T, typename... Args>
void raw_debug(istream_iterator<string> it, T a, Args... args) { cerr <<"<"<< *it << "->" << a << "> "; raw_debug(++it, args...); }
#define idebug(v) {cerr<<'['<<#v<<']';for(const auto &el:v)cerr << ' ' << el; cerr << endl;}
#define adebug(ar,n) {cerr<<'['<<#ar<<']';for(int my_imp_ind=0;my_imp_ind<n;++my_imp_ind)cerr << ' ' << ar[my_imp_ind]; cerr << endl;}
template <typename t1,typename t2> ostream &operator<<(ostream &os, const pair<t1,t2> &pr){return os<<"("<<pr.first<<";"<<pr.second<<")";};

#define rep(i,strt,end) for(int i = strt ; i !=int(end) ; (int(strt)<int(end))?++i:--i )
#define rall(vec) vec.rbegin(), vec.rend()
#define all(vec) vec.begin(), vec.end()
#define eb emplace_back
#define pb push_back
#define pob pop_back
#define pf push_front
#define pof pop_front

// template<typename mint>
// static constexpr uint64_t get_pr() {
//     uint64_t _mod = mint::mod();
//     using u64 = uint64_t;
//     u64 ds[32] = {};
//     int idx = 0;
//     u64 m = _mod - 1;
//     for (u64 i = 2; i * i <= m; ++i) {
//       if (m % i == 0) {
//         ds[idx++] = i;
//         while (m % i == 0) m /= i;
//       }
//     }
//     if (m != 1) ds[idx++] = m;

//     uint64_t _pr = 3;
//     while (1) {
//       int flg = 1;
//       for (int i = 0; i < idx; ++i) {
//         u64 a = _pr, b = (_mod - 1) / ds[i], r = 1;
//         while (b) {
//           if (b & 1) r = r * a % _mod;
//           a = a * a % _mod;
//           b >>= 1;
//         }
//         if (r == 1) {
//           flg = 0;
//           break;
//         }
//       }
//       if (flg == 1) break;
//       ++_pr;
//     }
//     return _pr;
//   };


mt19937_64 rng_64( chrono::steady_clock::now().time_since_epoch().count() );
constexpr int ilog2( int num ) { return 8*sizeof(int) - __builtin_clz( num ) - 1; }
template<typename tpow,typename texp=lint> constexpr tpow mpow(tpow x,unsigned long long e,tpow m){tpow res=1;while(e){if(e&1)res=(res*texp(1)*x)%m;e>>=1;x=(texp(x)*x)%m;}return res;}
// const unsigned long long my_mod=(27ull<<59)+1;
const unsigned long long my_mod=(549755813881ull<<24)+1;
const int mod=1e5;
#include"../modulo_int.cpp"
#include"../formalPowerSeries.cpp"
const int omod=998244353;
typedef modulo_int<my_mod,unsigned long long,__uint128_t> bmint;
typedef FormalPowerSeries<bmint> bfps;
template<typename tfps>
bool eq(const FormalPowerSeries<tfps> &A, const FormalPowerSeries<tfps> &B){
	const int n=max<int>(A.size(),B.size());
	bool res=1;rep(i,0,n){
		int va=A[i],vb=(int(B.size())>i)?int(B[i]):0;
		res&=va==vb;
	} return res;
}
fps my_mult(const fps&A,const fps&B){
	const int n=A.size()+B.size()-1,e=ilog2(n-1)+1;
	vector<bmint> a(1<<e,bmint(0ull)),b(1<<e,bmint(0ull));
	// vector<unsigned long long> vec(n,0);
	// rep(i,0,A.size())rep(j,0,B.size())vec[i]+=int(A[i])*1ll*int(B[i]);
	// rep(i,0,n)assert(vec[i]<bmint::mod());

	rep(i,0,A.size())a[i]=(unsigned long long)(A[i]);
	rep(i,0,B.size())b[i]=(unsigned long long)(B[i]);
	// idebug(a);idebug(b);
	bfps ab=bfps::mult_naive(bfps(a),bfps(b));
	bfps mab=bfps(a)*bfps(b);
	internal::fft<bmint>(a,0);internal::fft<bmint>(b,0);
	rep(i,0,1<<e)a[i]*=b[i];
	internal::fft<bmint>(a,1);
	// idebug(ab);
	// idebug(mab);
	// debug(eq<bmint>(ab,mab));
	// idebug(a);
	// rep(i,0,1<<e)assert(ab[i]==a[i]);
	// unsigned long long mx=0;rep(i,0,1<<e)mx=max<unsigned long long>(mx,(unsigned long long)a[i]);
	// debug("$#$$$$$$$$$$$$$$$$$$$$$",mx,bmint::mod());
	fps AB(n);rep(i,0,n)AB[i]=mint( (unsigned long long)(a[i]) );
	return AB;
} template<typename tint> vector<tint> gen_rand(int n){
	vector<tint> rs(n);rep(i,0,n)rs[i]=lint(rng_64());
	return rs;
} 
// static constexpr int level =  __builtin_ctzll(bmint::mod() - 1);
//   bmint dw[level], dy[level];

//   void setwy(int k) {
//   	int pr=3;
//     bmint w[level], y[level];
//     w[k - 1] = bmint(pr).pow((bmint::mod() - 1) / (1ull << k));
//     y[k - 1] = w[k - 1].inv();
//     for (int i = k - 2; i > 0; --i)
//       w[i] = w[i + 1] * w[i + 1], y[i] = y[i + 1] * y[i + 1];
//   	for(int i=0;i<k;++i)cerr << w[i] << " \n"[i+1==k];
//     dw[1] = w[1], dy[1] = y[1], dw[2] = w[2], dy[2] = y[2];
//     for (int i = 3; i < k; ++i) {
//       dw[i] = dw[i - 1] * y[i - 2] * w[i];
//       dy[i] = dy[i - 1] * w[i - 2] * y[i];
//     }
//   }

 void solve() {
 	// constexpr int r1=internal::countr_zero_constexpr(bmint::mod()-1);
    // std::array<bmint,r1+1> root; // precompute the roots
    // root[r1]=bmint(internal::primitive_root(bmint::mod())).pow((bmint::mod()-1)>>r1);
    // debug(get_pr<bmint>());setwy(level);
    // rep(i,0,level)cerr << dw[i] << " \n"[i+1==level];
    // for(int i=r1;i>0;--i)root[i-1]=root[i]*root[i];
	// rep(i,0,r1+1)cerr << root[i] << " \n"[i==r1];
	// rep(i,2,r1+1){
	// 	unsigned long long ut=mpow<unsigned long long,__uint128_t>(root[i],((unsigned long long)(1))<<(i-0),bmint::mod());
	// 	debug(ut,i,root[i]);
	// 	assert(ut==1ull);
	// }
 	// rep(i,0,r1+1)cerr << internal::rank_fft_root<bmint>[i] << " \n"[i==r1];
 	// debug(r1,bmint(-1));
 	// debug((bmint::mod()-1)>>r1,bmint(3).pow(27),internal::primitive_root(bmint::mod()));
 	// constexpr int r2=internal::countr_zero_constexpr(omod-1)-1;
 	// rep(i,0,r2+1)cerr << internal::rank_fft_root<modulo_int<omod>>[i] << " \n"[i==r2];
 	// rep(i,2,r2+1){
 	// 	modulo_int<omod> tmp=internal::rank_fft_root<modulo_int<omod>>[i];
	// 	int ut=mpow(int(tmp),1<<(i-0),omod);
	// 	debug(ut,i,tmp);
	// 	assert(ut==1);
	// }
 	// debug(r2,omod>>r2);
 	// return;

	const int n=1e3,iter=100;
	rep(i,0,iter){
		bfps A(gen_rand<bmint>(n)),B(gen_rand<bmint>(n));
		bfps bAB=bfps::mult_naive(A,B),AB=A*B;
		const bool bmul=eq(bAB,AB);

		fps C(gen_rand<mint>(n)),D(gen_rand<mint>(n));
		fps bCD=fps::mult_naive(C,D),mCD=my_mult(C,D);
		const bool mmul=eq(bCD,mCD);
		if(!bmul||!mmul){
			debug("error",bmul,mmul);if(!bmul){
				idebug(A);idebug(B);
				idebug(bAB);idebug(AB);
			} if(!mmul){
				idebug(C);idebug(D);
				idebug(bCD);idebug(mCD);
			} return;
		}
		debug(i);
	}debug("horaay");
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