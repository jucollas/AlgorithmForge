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
Does this all extend to transposed????
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

mt19937_64 rng_64( chrono::steady_clock::now().time_since_epoch().count() );
constexpr int ilog2( int num ) { return 8*sizeof(int) - __builtin_clz( num ) - 1; }
template<typename tpow,typename texp=lint> constexpr tpow mpow(tpow x,lint e,tpow m){tpow res=1;while(e){if(e&1ll)res=(res*texp(1)*x)%m;e>>=1;x=(x*texp(1)*x)%m;}return res;}
const int mod = 998244353;
#include"../../../modulo_int.cpp"
#include"fft.cpp"

/*
    vector<mint> B(2<<lgi,0);
    rep(i,0,n)A[i]*=x.pow(-(i*1ll*(i-1)/2ll));
    rep(i,0,2*n)B[i]=x.pow(i*1ll*(i-1)/2ll);
    reverse(A.begin(),A.end());A.resize(2<<lgi,0);
    internal::fft(A,0);internal::fft(B,0);
    rep(i,0,2<<lgi)A[i]*=B[i];
    internal::fft(A,1);A=vector<mint>(A.begin()+n-1,A.begin()+2*n-1);
    rep(i,0,n)A[i]*=x.pow(-(i*1ll*(i-1)/2ll));
    return A;
*/

vector<mint> brute_conv(const vector<mint>&v1,const vector<mint>&v2){
    const int n=v1.size();vector<mint> res(n,0);
    rep(i,0,n)rep(j,0,n)res[(i+j)%n]+=v1[i]*v2[j];
    return res;
}

vector<mint> fft_conv(vector<mint>F,vector<mint>G){
    // Inspiration from the toeplitz matrix part on 
    // "Computational frameworks for the fast fourier transform"
    const int n=F.size(),lgi=ilog2(n-1)+2; if( (n&(n-1))==0 ){
        internal::fft(F,0);internal::fft(G,0);
        rep(i,0,n)F[i]*=G[i];
        internal::fft(F,1);return F;
    } F.resize(1<<lgi,0);G.resize(1<<lgi,0);
    rep(i,0,n)F[(1<<lgi)-n+i]=F[i];
    internal::fft(F,0);internal::fft(G,0);
    rep(i,0,1<<lgi)F[i]*=G[i];
    internal::fft(F,1);
    F.resize(n); return F;
} vector<mint> chirpz(vector<mint>A,const mint x=internal::primitive_root(mod) ){
    // some stuff adapted from * https://nyaannyaan.github.io/library/ntt/chirp-z.hpp
    // and * https://codeforces.com/blog/entry/83532
    const int n=A.size(),lgi=ilog2(n-1)+1;
    vector<mint> B(2<<lgi,0),ixs(n);
    const mint ix=x.inv(); mint ixa=1,xa=1;
    ixs[0]=1; rep(i,1,n)ixs[i]=ixs[i-1]*ixa,ixa*=ix,A[i]*=ixs[i];
    B[0]=1; rep(i,1,2*n)B[i]=B[i-1]*xa,xa*=x;
    reverse(A.begin(),A.end());A.resize(2<<lgi,0);
    internal::fft(A,0);internal::fft(B,0);
    rep(i,0,2<<lgi)A[i]*=B[i];
    internal::fft(A,1);A=vector<mint>(A.begin()+n-1,A.begin()+2*n-1);
    rep(i,0,n)A[i]*=ixs[i];
    return A;
}
vector<mint> brute_chirpz(const vector<mint>&A,const mint x=internal::primitive_root(mod)){
    const int n=A.size();
    vector<mint> B(n); mint xx=1;rep(i,0,n){
        mint xxx=1;rep(j,0,n){
            B[i]+=A[j]*xxx;xxx*=xx;
        } xx*=x;
    } return B;
} inline int gen(int lim){return (lim+rng_64()%lim)%lim;}
inline vector<mint> vgen(int n,int lim){vector<mint> v(n);for(auto&ac:v)ac=gen(lim);return v;}
void solve() {
	const int n=3e3,iter=100,lim=mod;
    vector<mint>A,B;rep(i,0,iter){debug(i);
        A=vgen(n,lim);B=vgen(n,lim);
        vector<mint> rbz=brute_chirpz(A);
        vector<mint> rfz=chirpz(A);
        // idebug(A);idebug(rbz);idebug(rfz);
        // debug(rbz==rfz); 
        assert(rbz==rfz);

        vector<mint>rb=brute_conv(A,B);
        vector<mint>rf=fft_conv(A,B);
        
        // idebug(rb);idebug(rf);idebug(A);idebug(B);
        // debug(rb==rf);
        assert(rb==rf);
    }
    cerr << "temrino exitosamente" << endl;
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