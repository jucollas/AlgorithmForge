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
#define pb push_back
#define pob pop_back
#define pf push_front
#define pof pop_front

mt19937_64 rng_64( chrono::steady_clock::now().time_since_epoch().count() );
constexpr int ilog2( int num ) { return 8*sizeof(int) - __builtin_clz( num ) - 1; }
template<typename tpow> constexpr tpow mpow(tpow x,lint e,tpow m){tpow res=1;while(e){if(e&1ll)res=(res*1ll*x)%m;e>>=1;x=(x*1ll*x)%m;}return res;}


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

GCD_CONV tested in https://judge.yosupo.jp/problem/gcd_convolution
LCM_CONV tested in https://codeforces.com/gym/105053/problem/G
				   https://judge.yosupo.jp/problem/lcm_convolution
Note it wont work well for A*A

Taken from https://codeforces.com/blog/entry/119082

Im assuming from my template:
#define rep(i,strt,end) for(int i = strt ; i !=int(end) ; (int(strt)<int(end))?++i:--i )
*/


template<typename T>
class div_conv{
	vector<int> prm;
public:
	constexpr div_conv(int limit){
		vector<bool> crb(limit,1);crb[0]=crb[1]=0;
		rep(i,0,limit)if(crb[i]){
			prm.push_back(i);
			if(i*1ll*i<lint(limit))for(int j=i*i;j<limit;j+=i){
				crb[j]=0;
			}
		}
	}

	void zeta_mult(vector<T> &vc)const{
		// $vc'_x=\sum_{x|y}vc_y$
		const int n=(vc.size()-1);
		for(int p:prm){
			if(p>n)break;
			for(int i=n/p;i;--i) vc[i]+=vc[i*p];
		}
	}
	void mobi_mult(vector<T> &vc)const{
		// $vc_x=\sum_{x|y}vc'_y$
		const int n=(vc.size()-1);
		for(int p:prm){
			if(p>n)break;
			for(int i=1;i*p<=n;++i)vc[i]-=vc[i*p];
		}
	}
	void gcd_conv(vector<T> &A,vector<T> &B)const{
		// I assume |A|=|B|
		// $A'_x=\sum_{x=gcd(u,v)}A_uB_v$
		zeta_mult(A);zeta_mult(B);
		rep(i,0,A.size())A[i]*=B[i];
		mobi_mult(A);
	}
	
	void zeta_div(vector<T> &vc)const{
		// $vc_x=\sum_{y|x}vc'_y$
		const int n=(vc.size()-1);
		for(int p:prm){
			if(p>n)break;
			for(int i=1;i*p<=n;++i)vc[i*p]+=vc[i];
		}
	}
	void mobi_div(vector<T> &vc)const{
		// $vc'_x=\sum_{y|x}vc_y$
		const int n=(vc.size()-1);
		for(int p:prm){
			if(p>n)break;
			for(int i=n/p;i;--i) vc[i*p]-=vc[i];
		}
	}
	
	void lcm_conv(vector<T> &A,vector<T> &B)const{
		// I assume |A|=|B|
		// $A'_x=\sum_{x=lcm(u,v)}A_uB_v$
		zeta_div(A);zeta_div(B);
		rep(i,0,A.size())A[i]*=B[i];
		mobi_div(A);
	}

};
const int limit =1e6;
const div_conv<mint> dc(limit);

void solve() {
	auto gen=[&](){int tmp = rng_64()%(limit-2);if(tmp<0)tmp+=limit-2;++tmp;return tmp;};
	const int n=1e4; static_assert(n<=limit);
	vector<mint> as,bs;
	vector<mint> A(n); rep(i,1,n)A[i]=rng_64();
	vector<mint> B(n); rep(i,1,n)B[i]=rng_64();
	vector<mint> C(n,0);
	rep(i,1,n)rep(j,1,n)C[__gcd(i,j)]+=A[i]*B[j];
	as=A;bs=B;
	dc.gcd_conv(A,B);
	if(A!=C){
		idebug(C);idebug(A);
		idebug(as);idebug(bs);
	}else cout << "ok"<<endl;
}

int32_t main(){
	ios_base::sync_with_stdio(false);
    cin.tie(NULL);
	cout << setprecision(12) << fixed;

    int t = 10;
    // cin >> t; ++t;
    while ( --t ) {
		solve();
    }
	return 0;
}

