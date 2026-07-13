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

#include <bits/stdc++.h>

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
template<typename tpow> constexpr tpow mpow(tpow x,lint e,tpow m){tpow res=1;while(e){if(e&1ll)res=(res*1ll*x)%m;e>>=1;x=(x*1ll*x)%m;}return res;}

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

template<int m>
struct modulo_int{ static_assert(m>0);
	constexpr static int mod(){return m;}
	
	int vl;
	constexpr modulo_int()noexcept:vl(0){};
	constexpr modulo_int( int v)noexcept:vl(v>=0?(v<m?v:v%m):(v+m>=0?v+m:(v%m)+m)){};
	constexpr modulo_int(lint v)noexcept:vl(v>=0?(v<m?v:v%m):(v+m>=0?v+m:(v%m)+m)){};
	constexpr modulo_int(unsigned long long v)noexcept:vl(v<m?v:v%m){};
	
	modulo_int &operator +=(const modulo_int &ot){ vl+=ot.vl; if(vl>=m)vl-=m; return *this; }
	modulo_int  operator + (const modulo_int &ot)const{ return modulo_int(*this)+=ot; }
	modulo_int &operator -=(const modulo_int &ot){ vl-=ot.vl; if(vl<0)vl+=m; return *this; }
	modulo_int  operator - (const modulo_int &ot)const{ return modulo_int(*this)-=ot; }
	modulo_int &operator *=(const modulo_int &ot){ vl=(vl*1ll*ot.vl)%m; return *this; }
	modulo_int  operator * (const modulo_int &ot)const{ return modulo_int(*this)*=ot; }
	modulo_int &operator /=(const modulo_int &ot){ (*this)*=ot.inv(); return *this; }
	modulo_int  operator / (const modulo_int &ot)const{ return modulo_int(*this)/=ot; }
	
	modulo_int inv()const{return modulo_int(vl).pow(m-2);}//Fermats little theorem
	modulo_int operator -()const {return modulo_int(-vl);}
	modulo_int pow(lint e)const{return modulo_int(mpow(vl,e%(m-1),m));}
	
	bool operator ==(const modulo_int &ot)const{return vl==ot.vl;} 
	bool operator ==(const  int &ot)const{return vl==ot;   }
	bool operator !=(const modulo_int &ot)const{return vl!=ot.vl;}
	
	operator bool() const { return vl; }
	operator  int() const { return vl; }
	
	friend ostream &operator<<(ostream &os,const modulo_int &ac){return os << ac.vl;}
	friend istream &operator>>(istream&is,modulo_int &ac){int v;is>>v;ac=modulo_int(v);return is;}	
}; typedef modulo_int<mod> mint;
/*
Autor: Oscar Vargas Pabon

and convolution tested in https://judge.yosupo.jp/problem/bitwise_and_convolution
or convolution currently untested

Note it wont work well for A*A

Taken from https://codeforces.com/blog/entry/119082
			https://codeforces.com/blog/entry/92128
			https://gist.github.com/dario2994/e3257326ee80c054d3b48766b600991a

Im assuming from my template:
#define rep(i,strt,end) for(int i = strt ; i !=int(end) ; (int(strt)<int(end))?++i:--i )
*/

template<typename T>
void trans_subset(vector<T> &vec,int sd){
	// I assume $sd\in\{-1,1\}$
	// When sd=+1 $vec[k]'=\sum_{i,j\subseteq k}vec[i] *vec[j] $
	// When sd=-1 $vec[k] =\sum_{i,j\subseteq k}vec[i]'*vec[j]'$
	for(int e=1;e<int(vec.size());e*=2)rep(i,0,vec.size()){
		if(i&e)vec[i]+=vec[i^e]*T(sd);
	}
}
template<typename T>
vector<T> or_conv(vector<T> A,vector<T> B){
	// the answer is shown in A; I assume $|A|=|B|=2^x$ for some x
	// Computes $A'[k]=\sum_{(i|j)==k}A[i]*B[j]$
	trans_subset(A,1);trans_subset(B,1);
	rep(i,0,A.size())A[i]*=B[i];
	trans_subset(A,-1);
	return A;
}
template<typename T>	
void trans_superset(vector<T> &vec,int sd){
	// I assume $sd\in\{-1,1\}$
	// When sd=+1 $vec[k]'=\sum_{i,j\superseteq k}vec[i] *vec[j] $
	// When sd=-1 $vec[k] =\sum_{i,j\superseteq k}vec[i]'*vec[j]'$
	for(int e=1;e<int(vec.size());e*=2)rep(i,0,vec.size()){
		if(i&e)vec[i^e]+=vec[i]*T(sd);
	}
}
template<typename T>	
vector<T> and_conv(vector<T> &A,vector<T> &B){
	// the answer is shown in A; I assume $|A|=|B|=2^x$ for some x
	// Computes A'[k]=\sum_{(i&j)==k}A[i]*B[j]$
	trans_superset(A,1);trans_superset(B,1);
	rep(i,0,A.size())A[i]*=B[i];
	trans_superset(A,-1);
	return A;
}
template<typename T>
vector<T> subset_conv(const vector<T> &A,const vector<T>&B){
	// the answer is shown in A; I assume $|A|=|B|=2^x$ for some x
	// Computes $C[k]=\sum_{s\subseteq k}A[s]*B[k\setminus s]$
	const int n=ilog2(A.size()); assert(A.size()==B.size()&&(1<<n)==int(A.size()));
	vector<vector<T>> A_hat(n+1,vector<T>(1<<n,0)),B_hat=A_hat;
	rep(i,0,1<<n){ int pcnt=__builtin_popcount(i);
		A_hat[pcnt][i]=A[i]; B_hat[pcnt][i]=B[i];
	} rep(i,0,n+1)trans_subset(A_hat[i],1),trans_subset(B_hat[i],1);
	vector<T> C(1<<n),C_hat(1<<n);rep(k,0,n+1){
		fill(all(C_hat),T(0));rep(i,0,k+1)rep(j,0,1<<n){
			C_hat[j]+=A_hat[i][j]*B_hat[k-i][j];
		}trans_subset(C_hat,-1);
		rep(i,0,1<<n)if(__builtin_popcount(i)==k)C[i]=C_hat[i];
	}return C;
}
template<typename T> vector<T> subset_iconv(const vector<T>&A,const vector<T>&C){
	// Computes B such that $C[k]=\sum_{i\subseteq k}A[i]*B[k\setminus i]$
	const T A0_inv=A[0].inv(); assert(!(A[0]==0));// Im assuming some .inv()
	const int n=ilog2(A.size()); assert(A.size()==C.size()&&(1<<n)==int(A.size()));
	vector<vector<T>> A_hat(n+1,vector<T>(1<<n,0)),B_hat=A_hat;
	rep(i,0,1<<n)A_hat[__builtin_popcount(i)][i]=A[i];
	rep(i,0,n+1)trans_subset(A_hat[i],1);
	vector<T>B(1<<n); rep(k,0,n+1){
		rep(i,0,k)rep(j,0,1<<n)B_hat[k][j]+=B_hat[i][j]*A_hat[k-i][j];
		trans_subset(B_hat[k],-1);
		rep(j,0,1<<n){
			// C_hat[k][j] = B_hat[k] =0 since ppcnt(j)!=k
			if(__builtin_popcount(j)!=k)B_hat[k][j]=0;
			//C_hat[k][j]=C[j]=B[j]*A[0]+B_hat[k][j]
			else B[j]=B_hat[k][j]=(C[j]-B_hat[k][j])*A0_inv;
		} trans_subset(B_hat[k],1);
	} return B;
}
template<typename T> vector<T> set_exp(const vector<T>&A){
	// Computes $exp(A)=\sum_{i\geq 0}A^i/i!$
	// exp(A)[k]=\sum_{p\in P[k]}\prod_{i\in p}A_{p_i}$
	// Where P is the set of all partitions of k into nonempty sets
	assert(A[0]==0);// required for A^i to be nilpotent
	vector<T> B={1};B.reserve(A.size());
	for(int e=1;e<int(A.size());e<<=1){
		for(const T &ac:subset_conv({A.begin()+e,A.begin()+e*2},B)){
			B.emplace_back(ac);
		}
	}return B;//B=exp(A)
}
template<typename T> vector<T> set_log(const vector<T>&A){
	// Computes given A=exp(B) ; computes B;
	assert(A[0]==1);
	vector<T>B={0};B.reserve(A.size());
	for(int e=1;e<int(A.size());e<<=1){
		for(const T &ac:subset_iconv<T>( {A.begin(),A.begin()+e},{A.begin()+e,A.begin()+e*2} )){
			B.emplace_back(ac);
		}
	}return B;
}


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

