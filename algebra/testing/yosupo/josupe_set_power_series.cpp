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


this test has 
	https://judge.yosupo.jp/submission/349368
The set_power_series impl has
	https://judge.yosupo.jp/submission/349366
*/
#include <cassert>
#include <bits/stdc++.h>
#pragma GCC optimize("O3")
#define NDEBUG

typedef long long lint;

using namespace std;

#define debug(args...) { string _s = #args; replace(_s.begin(), _s.end(), ',', ' '); stringstream _ss(_s); istream_iterator<string> _it(_ss); raw_debug(_it, args);}
void raw_debug(istream_iterator<string> it) {cerr<<endl;assert(it==it);}
template<typename T, typename... Args>
void raw_debug(istream_iterator<string> it, T a, Args... args) { cerr <<"<"<< *it << "->" << a << "> "; raw_debug(++it, args...); }
#define idebug(v) {cerr<<'['<<#v<<']';for(const auto &el:v)cerr << ' ' << el; cerr << endl;}
#define adebug(ar,n) {cerr<<'['<<#ar<<']';for(int i=0;i<n;++i)cerr << ' ' << ar[i]; cerr << endl;}

#define rep(i,strt,end) for(int i = strt ; i !=int(end) ; (int(strt)<int(end))?++i:--i )
#define rall(vec) vec.rbegin(), vec.rend()
#define all(vec) vec.begin(), vec.end()
#define sz(vec) int(vec.size())
#define eb emplace_back
#define pb push_back
#define pob pop_back
#define pf push_front
#define pof pop_front

mt19937_64 rng_64( chrono::steady_clock::now().time_since_epoch().count() );
constexpr int ilog2( int num ) { return 8*sizeof(int) - __builtin_clz( num ) - 1; }
template<typename tpow,typename texp=int64_t> constexpr tpow mpow(tpow x,uint64_t e,tpow m){tpow res=1;while(e){if(e&1)res=(texp(res)*x)%m;e>>=1;x=(texp(x)*x)%m;}return res;}

// const int template_limit = 1e6;
// int a[template_limit], b[template_limit];
const int mod = 998244353;
/*
Author: Oscar Vargas Pabon

It is though to work modulo primes, so inv (.inv,/,/=) and pow may not
	work properly otherwise.

I assume from my template :
inv :: int mpow(int x,int e,int m){int res=1;while(e){if(e&1)res=(res*1ll*x)%m;e>>=1;x=(x*1ll*x)%m;}return res;}	

Tested in testing/test_alghelp.cpp and in fft/ntt stuff
*/
template<uint64_t raw_m,typename tint=uint32_t,typename tmul=uint64_t,bool arbi_ntt=0>
struct modulo_int{ constexpr static tint m=raw_m; static_assert(m>0);
	constexpr static tint mod(){return m;}
	constexpr static bool arbitrary_ntt(){return arbi_ntt;}
	
	tint vl;
	inline constexpr modulo_int()noexcept:vl(0){};
	inline constexpr modulo_int(      int v)noexcept:vl(v>=0?(v<m?v:v%m):(v+m>=0?v+m:(m+(v%m))%m)){};
	inline constexpr modulo_int(long long v)noexcept:vl(v>=0?(v<m?v:v%m):(v+m>=0?v+m:(m+(v%m))%m)){};
	inline constexpr modulo_int(uint32_t v)noexcept:vl(v<m?v:v%m){};
	inline constexpr modulo_int(uint64_t v)noexcept:vl(v<m?v:v%m){};
	
	inline constexpr modulo_int &operator +=(const modulo_int &ot){ vl= m-vl>ot.vl?vl+ot.vl:ot.vl-(m-vl); return *this; }
	inline constexpr modulo_int  operator + (const modulo_int &ot)const{ return modulo_int(*this)+=ot; }
	inline constexpr modulo_int &operator -=(const modulo_int &ot){ vl=(vl>=ot.vl)?vl-ot.vl:(m-ot.vl)+vl; return *this; }
	inline constexpr modulo_int  operator - (const modulo_int &ot)const{ return modulo_int(*this)-=ot; }
	inline constexpr modulo_int &operator *=(const modulo_int &ot){ vl=(tmul(vl)*ot.vl)%m; return *this; }
	inline constexpr modulo_int  operator * (const modulo_int &ot)const{ return modulo_int(*this)*=ot; }
	inline constexpr modulo_int &operator /=(const modulo_int &ot){ (*this)*=ot.inv(); return *this; }
	inline constexpr modulo_int  operator / (const modulo_int &ot)const{ return modulo_int(*this)*=ot.inv(); }
	
	inline constexpr modulo_int operator -()const {return vl?m-vl:0;}
	inline constexpr modulo_int inv()const{return mpow<tint,tmul>(vl,m-2,m);}//Fermats little theorem
	inline constexpr modulo_int pow(long long e)const{return e>=0?modulo_int(mpow<tint,tmul>(vl,e%(m-1),m)):inv().pow(-e);}
	
	inline constexpr bool operator ==(const modulo_int &ot)const{return vl==ot.vl;} 
	inline constexpr bool operator ==(const  int &ot)const{return vl==ot;}
	inline constexpr bool operator !=(const modulo_int &ot)const{return vl!=ot.vl;}
	
	inline constexpr operator     bool() const{return vl;}
	inline constexpr operator  int32_t() const{return vl;}
	inline constexpr operator uint32_t() const{return vl;}
	inline constexpr operator  int64_t() const{return vl;}
	inline constexpr operator uint64_t() const{return vl;}
	
	friend ostream &operator<<(ostream &os,const modulo_int &ac){return os << ac.vl;}
	friend istream &operator>>(istream&is,modulo_int &ac){int v;is>>v;ac=modulo_int(v);return is;}	
}; typedef modulo_int<mod> mint;

/*
Author: Oscar Vargas Pabon

and convolution tested in https://judge.yosupo.jp/problem/bitwise_and_convolution
or convolution currently untested

exp and log tested in https://judge.yosupo.jp/submission/349366

All impl are also tested in my personal libs

Taken from https://codeforces.com/blog/entry/119082
			https://codeforces.com/blog/entry/92128
			https://gist.github.com/dario2994/e3257326ee80c054d3b48766b600991a

Im assuming from my template:
#define rep(i,strt,end) for(int i = strt ; i !=int(end) ; (int(strt)<int(end))?++i:--i )
*/
/*
	for(int e=1;e<n;e*=2)for(int i=0;i<n;++i)if(i&e){
		if constexpr(tp)vec[i^e]+=vec[i];
		else            vec[i^e]-=vec[i]; }
*/
const int max_exp=22;// (1<<max_exp)> maximum array size possible
// #define VER1
template<typename tint>void raw_subset(tint*vec,int n,bool tp){
#ifdef VER1
	if(tp){
	for(int e=1;e<n;e*=2)for(int i=0;i<n;++i)if(i&e)vec[i]+=vec[i^e];
	} else {
	for(int e=1;e<n;e*=2)for(int i=0;i<n;++i)if(i&e)vec[i]-=vec[i^e];
	}
#else
	if(tp){
	for(int e=1;e<n;e*=2)for(int r=0;r<n;r+=2*e)for(int l=0;l<e;++l)
		vec[l|r|e]+=vec[l|r];
	} else {
	for(int e=1;e<n;e*=2)for(int r=0;r<n;r+=2*e)for(int l=0;l<e;++l)
		vec[l|r|e]-=vec[l|r];
	}
#endif
}template<typename tint>inline void subset(vector<tint>&vec,bool tp){
	// When tp=1 $vec[k]'=\sum_{i,j\subseteq k}vec[i] *vec[j] $
	// When tp=0 $vec[k] =\sum_{i,j\subseteq k}vec[i]'*vec[j]'$
	raw_subset<tint>(vec.data(),vec.size(),tp);
}template<typename tint>void raw_superset(tint*vec,int n,bool tp){
#ifdef VER1
	if(tp){
	for(int e=1;e<n;e*=2)for(int i=0;i<n;++i)if(i&e)vec[i^e]+=vec[i];
	} else {
	for(int e=1;e<n;e*=2)for(int i=0;i<n;++i)if(i&e)vec[i^e]-=vec[i];
	}
#else
	if(tp){
	for(int e=1;e<n;e*=2)for(int r=0;r<n;r+=2*e)for(int l=0;l<e;++l)
		vec[l|r]+=vec[l|r|e];
	} else {
	for(int e=1;e<n;e*=2)for(int r=0;r<n;r+=2*e)for(int l=0;l<e;++l)
		vec[l|r]-=vec[l|r|e];
	}
#endif
}template<typename tint>inline void superset(vector<tint>&vec,bool tp){
	// When tp=1 $vec[k]'=\sum_{i,j\superseteq k}vec[i] *vec[j] $
	// When tp=0 $vec[k] =\sum_{i,j\superseteq k}vec[i]'*vec[j]'$
	raw_superset<tint>(vec.data(),vec.size(),tp);
} template<typename tint>vector<tint>or_conv(vector<tint>A,vector<tint>B){
	// the answer is shown in A; I assume $|A|=|B|=2^x$ for some x
	// Computes $A'[k]=\sum_{(i|j)==k}A[i]*B[j]$
	subset<tint>(A,1);subset<tint>(B,1);
	for(int i=0;i<int(A.size());++i)A[i]*=B[i];
	subset<tint>(A,0); return A;
} template<typename tint>vector<tint>and_conv(vector<tint>A,vector<tint>B){
	// the answer is shown in A; I assume $|A|=|B|=2^x$ for some x
	// Computes A'[k]=\sum_{(i&j)==k}A[i]*B[j]$
	superset<tint>(A,1);superset<tint>(B,1);
	for(int i=0;i<int(A.size());++i)A[i]*=B[i];
	superset<tint>(A,0); return A;
} template<typename tint>void raw_subset_conv(const tint A[],const tint B[],tint*C,int n){
	// I assume $2^n=|A|=|B|=|C|$
	static tint A_hat[(max_exp+1)<<max_exp],B_hat[(max_exp+1)<<max_exp],C_hat[(max_exp+1)<<max_exp];
	for(int i=0;i<((n+1)<<n);++i)A_hat[i]=B_hat[i]=C_hat[i]=0;
	for(int i=0;i<(1<<n);++i){int pcnt=__builtin_popcount(i);
		A_hat[pcnt<<n|i]=A[i]; B_hat[pcnt<<n|i]=B[i];
	}for(int i=0;i<=n;++i)raw_subset<tint>(A_hat+(i<<n),1<<n,1),
	                      raw_subset<tint>(B_hat+(i<<n),1<<n,1);
	for(int k=0;k<=n;++k)for(int i=0;i<=k;++i)for(int j=0;j<(1<<n);++j)
		C_hat[k<<n|j]+=A_hat[i<<n|j]*B_hat[(k-i)<<n|j];
	for(int i=0;i<=n;++i)raw_subset<tint>(C_hat+(i<<n),1<<n,0);
	for(int i=0;i<(1<<n);++i)C[i]=C_hat[__builtin_popcount(i)<<n|i];
}template<typename tint>inline vector<tint>subset_conv(const vector<tint>&A,const vector<tint>&B){
	// Computes $C[k]=\sum_{s\subseteq k}A[s]*B[k\setminus s]$
	const int n=ilog2(A.size()); assert(A.size()==B.size()&&(1<<n)==int(A.size()));
	vector<tint> C(1<<n); raw_subset_conv<tint>(A.data(),B.data(),C.data(),n);
	return C;
} template<typename tint>void raw_subset_iconv(const tint A[],const tint C[],tint B[],int n){
	// I assume $2^n=|A|=|B|=|C|$
	const tint A0_inv=A[0].inv(); assert(!(A[0]==0));// Im assuming some .inv()
	static tint A_hat[(max_exp+1)<<max_exp],B_hat[(max_exp+1)<<max_exp];
	for(int i=0;i<((n+1)<<n);++i)A_hat[i]=B_hat[i]=0;
	for(int i=0;i<(1<<n);++i)A_hat[__builtin_popcount(i)<<n|i]=A[i];
	for(int i=0;i<=n;++i)raw_subset<tint>(A_hat+(i<<n),1<<n,1);
	for(int k=0;k<=n;++k){
		for(int i=0;i<k;++i)for(int j=0;j<(1<<n);++j)
			B_hat[k<<n|j]+=B_hat[i<<n|j]*A_hat[(k-i)<<n|j];
		raw_subset<tint>(B_hat+(k<<n),1<<n,0);
		for(int j=0;j<(1<<n);++j){
			if(__builtin_popcount(j)!=k)B_hat[k<<n|j]=0;
			else B[j]=B_hat[k<<n|j]=(C[j]-B_hat[k<<n|j])*A0_inv;
		} raw_subset<tint>(B_hat+(k<<n),1<<n,1); }
}template<typename tint>inline vector<tint>subset_iconv(const vector<tint>&A,const vector<tint>&C){
	// Computes B such that $C[k]=\sum_{i\subseteq k}A[i]*B[k\setminus i]$
	const int n=ilog2(A.size()); assert(A.size()==C.size()&&(1<<n)==int(A.size()));
	vector<tint>B(1<<n);raw_subset_iconv<tint>(A.data(),C.data(),B.data(),n);return B;
} template<typename tint>vector<tint>set_exp(const vector<tint>&A){
	// Computes $exp(A)=\sum_{i\geq 0}A^i/i!$
	// exp(A)[k]=\sum_{p\in P[k]}\prod_{i\in p}A_{p_i}$
	// Where P is the set of all partitions of k into nonempty sets
	assert(A[0]==0);// required for A^i to be nilpotent
	vector<tint> B(A.size());B[0]=1;
	for(int e=0;(1<<e)<int(A.size());++e)
		raw_subset_conv<tint>(A.data()+(1<<e),B.data(),B.data()+(1<<e),e);
	return B;//B=exp(A)
}template<typename tint>vector<tint>set_log(const vector<tint>&A){
	// Computes given A=exp(B) ; computes B;
	assert(A[0]==1);
	vector<tint>B(A.size(),0);
	for(int e=0;(1<<e)<int(A.size());++e)
		raw_subset_iconv<tint>(A.data(),A.data()+(1<<e),B.data()+(1<<e),e);
	return B; }
/* end set-stuff */

void solve() {
	int n;cin>>n;
	vector<mint> A(1<<n);rep(i,0,1<<n)cin>>A[i];
	vector<mint> B=set_exp<mint>(A);
	vector<mint> C=set_log<mint>(B);assert(C==A);
	for(mint ac:B)cout << ac << ' ';
	cout << '\n';
}

int32_t main(){
	ios_base::sync_with_stdio(false);
    cin.tie(NULL);
	cout << setprecision(12) << fixed;

    int t = 2;
    // cin >> t; ++t;
    while ( --t ) {
		solve();
    }
	return 0;
}

