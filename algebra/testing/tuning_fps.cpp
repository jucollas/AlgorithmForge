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
constexpr int mpow(int x,int e,int m){int res=1;while(e){if(e&1)res=(res*1ll*x)%m;e>>=1;x=(x*1ll*x)%m;}return res;}

const int template_limit = 1e6;
// int a[template_limit], b[template_limit];

const int mod = 998244353;

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
Author: Oscar Vargas Pabon

NTT taken from the atcoder library
	https://github.com/atcoder/ac-library/

Based on code by MarcosK and other people
	https://codeforces.com/contest/438/submission/340901913
and multiple blogs all around the place
	https://cp-algorithms.com/algebra/polynomial.html#inverse-series_1
	https://codeforces.com/blog/entry/56422
	https://codeforces.com/blog/entry/12513 - problem E

Tested in
	https://judge.yosupo.jp/problem/inv_of_formal_power_series
	https://judge.yosupo.jp/problem/exp_of_formal_power_series
	https://judge.yosupo.jp/problem/log_of_formal_power_series
	https://judge.yosupo.jp/problem/pow_of_formal_power_series
	https://judge.yosupo.jp/problem/sqrt_of_formal_power_series

The first version of this impl is in atcoder fps_24 A

I assume mint from algebra/mint.cpp

Dependencies I assume
everything     :: #define rep(i,strt,end) for(int i = strt ; i !=int(end) ; (int(strt)<int(end))?++i:--i )
everything     :: #define sz(vec) int(vec.size())
mult,square    :: int ilog2( int num ) { return 8*sizeof(int) - __builtin_clz( num ) - 1; }
modInverse     :: int mpow(int x,int e,int m){int res=1;while(e){if(e&1)res=(res*1ll*x)%m;e>>=1;x=(x*1ll*x)%m;}return res;}	
TonelliShanks  :: mt19937_64 rng_64( chrono::steady_clock::now().time_since_epoch().count() );
sqrt           :: TonelliShanks, tfps(?).inv()
pow            :: tfps(?).inv(), tfps(?).pow()



*/

/* START OF NTT */
namespace internal { // taken from atcoder

int countr_zero(unsigned int n) { return __builtin_ctz(n); }
constexpr int countr_zero_constexpr(unsigned int n) { int x = 0; while (!(n & (1 << x))) x++; return x; }

struct fft_info {
	static const int g =3; // primitive root 
	
    static constexpr int rank2 = countr_zero_constexpr(mint::mod() - 1);
    std::array<mint, rank2 + 1> root;   // root[i]^(2^i) == 1
    std::array<mint, rank2 + 1> iroot;  // root[i] * iroot[i] == 1

    std::array<mint, std::max(0, rank2 - 2 + 1)> rate2;
    std::array<mint, std::max(0, rank2 - 2 + 1)> irate2;

    std::array<mint, std::max(0, rank2 - 3 + 1)> rate3;
    std::array<mint, std::max(0, rank2 - 3 + 1)> irate3;

    fft_info() {
        root[rank2] = mint(g).pow((mint::mod() - 1) >> rank2);
        iroot[rank2] = root[rank2].inv();
        for (int i = rank2 - 1; i >= 0; i--) {
            root[i] = root[i + 1] * root[i + 1];
            iroot[i] = iroot[i + 1] * iroot[i + 1];
        }

        {
            mint prod = 1, iprod = 1;
            for (int i = 0; i <= rank2 - 2; i++) {
                rate2[i] = root[i + 2] * prod;
                irate2[i] = iroot[i + 2] * iprod;
                prod *= iroot[i + 2];
                iprod *= root[i + 2];
            }
        }
        {
            mint prod = 1, iprod = 1;
            for (int i = 0; i <= rank2 - 3; i++) {
                rate3[i] = root[i + 3] * prod;
                irate3[i] = iroot[i + 3] * iprod;
                prod *= iroot[i + 3];
                iprod *= root[i + 3];
            }
        }
    }
};


void butterfly(std::vector<mint>& a) {
    int n = int(a.size());
    int h = countr_zero((unsigned int)n);

    static const fft_info info;

    int len = 0;  // a[i, i+(n>>len), i+2*(n>>len), ..] is transformed
    while (len < h) {
        if (h - len == 1) {
            int p = 1 << (h - len - 1);
            mint rot = 1;
            for (int s = 0; s < (1 << len); s++) {
                int offset = s << (h - len);
                for (int i = 0; i < p; i++) {
                    auto l = a[i + offset];
                    auto r = a[i + offset + p] * rot;
                    a[i + offset] = l + r;
                    a[i + offset + p] = l - r;
                }
                if (s + 1 != (1 << len))
                    rot *= info.rate2[countr_zero(~(unsigned int)(s))];
            }
            len++;
        } else {
            // 4-base
            int p = 1 << (h - len - 2);
            mint rot = 1, imag = info.root[2];
            for (int s = 0; s < (1 << len); s++) {
                mint rot2 = rot * rot;
                mint rot3 = rot2 * rot;
                int offset = s << (h - len);
                for (int i = 0; i < p; i++) {
                    auto mod2 = 1ULL * mint::mod() * mint::mod();
                    auto a0 = 1ULL * a[i + offset].vl;
                    auto a1 = 1ULL * a[i + offset + p].vl * rot.vl;
                    auto a2 = 1ULL * a[i + offset + 2 * p].vl * rot2.vl;
                    auto a3 = 1ULL * a[i + offset + 3 * p].vl * rot3.vl;
                    auto a1na3imag =
                        1ULL * mint(a1 + mod2 - a3).vl * imag.vl;
                    auto na2 = mod2 - a2;
                    a[i + offset] = a0 + a2 + a1 + a3;
                    a[i + offset + 1 * p] = a0 + a2 + (2 * mod2 - (a1 + a3));
                    a[i + offset + 2 * p] = a0 + na2 + a1na3imag;
                    a[i + offset + 3 * p] = a0 + na2 + (mod2 - a1na3imag);
                }
                if (s + 1 != (1 << len))
                    rot *= info.rate3[countr_zero(~(unsigned int)(s))];
            }
            len += 2;
        }
    }
}

void butterfly_inv(std::vector<mint>& a) {
    int n = int(a.size());
    int h = countr_zero((unsigned int)n);

    static const fft_info info;

    int len = h;  // a[i, i+(n>>len), i+2*(n>>len), ..] is transformed
    while (len) {
        if (len == 1) {
            int p = 1 << (h - len);
            mint irot = 1;
            for (int s = 0; s < (1 << (len - 1)); s++) {
                int offset = s << (h - len + 1);
                for (int i = 0; i < p; i++) {
                    auto l = a[i + offset];
                    auto r = a[i + offset + p];
                    a[i + offset] = l + r;
                    a[i + offset + p] =
                        (unsigned long long)((unsigned int)(l.vl - r.vl) + mint::mod()) *
                        irot.vl;
                    ;
                }
                if (s + 1 != (1 << (len - 1)))
                    irot *= info.irate2[countr_zero(~(unsigned int)(s))];
            }
            len--;
        } else {
            // 4-base
            int p = 1 << (h - len);
            mint irot = 1, iimag = info.iroot[2];
            for (int s = 0; s < (1 << (len - 2)); s++) {
                mint irot2 = irot * irot;
                mint irot3 = irot2 * irot;
                int offset = s << (h - len + 2);
                for (int i = 0; i < p; i++) {
                    auto a0 = 1ULL * a[i + offset + 0 * p].vl;
                    auto a1 = 1ULL * a[i + offset + 1 * p].vl;
                    auto a2 = 1ULL * a[i + offset + 2 * p].vl;
                    auto a3 = 1ULL * a[i + offset + 3 * p].vl;

                    auto a2na3iimag =
                        1ULL *
                        mint((mint::mod() + a2 - a3) * iimag.vl).vl;

                    a[i + offset] = a0 + a1 + a2 + a3;
                    a[i + offset + 1 * p] =
                        (a0 + (mint::mod() - a1) + a2na3iimag) * irot.vl;
                    a[i + offset + 2 * p] =
                        (a0 + a1 + (mint::mod() - a2) + (mint::mod() - a3)) *
                        irot2.vl;
                    a[i + offset + 3 * p] =
                        (a0 + (mint::mod() - a1) + (mint::mod() - a2na3iimag)) *
                        irot3.vl;
                }
                if (s + 1 != (1 << (len - 2)))
                    irot *= info.irate3[countr_zero(~(unsigned int)(s))];
            }
            len -= 2;
        }
    }
}

void fft(std::vector<mint> &A,bool invert){
	if(invert){
		internal::butterfly_inv(A);
		mint n_1 = mint(sz(A)).inv();
		for (mint & x : A)x*=n_1;
	} else internal::butterfly(A);
}

} //end internal namespace

int Tonelli_Shanks(mint a) {
	// usado por sqrt cuando trabajo con enteros modulo algo
	//plagiado epicamente de https://judge.yosupo.jp/submission/270105
	
	if (a < 2) return a;
	if (mpow(a, (mod - 1) / 2, mod) != 1) return -1;
	if (mod % 4 == 3) return mpow(a, (mod + 1) / 4, mod);

	mint b = 3;
	if (mod != 998244353) {
		while (mpow(b, (mod - 1) / 2, mod) == 1) {
			b=mint(int(rng_64()%(mod-3)) + 2);
		}
	}

	int q = mod - 1,Q = 0;
	while ( !(q&1) ) Q++, q /= 2;

	mint x = mpow(a, (q + 1) / 2, mod);
	b = mpow(b, q, mod);

	int shift = 2;
	while ( x*x != a) {
		mint error= mint(mpow(a,mod-2,mod))*x*x;
		if (mpow(error, 1 << (Q - shift), mod) != 1) x *= b;
		b *=b;
		++shift;
	}
	return x;
}

template<typename tfps>
struct FormalPowerSeries{
	std::vector<tfps> F;
	static FormalPowerSeries<tfps> mult_naive(const FormalPowerSeries<tfps>&A,
												const FormalPowerSeries<tfps>&B){
		FormalPowerSeries<tfps> C(sz(A)+sz(B),0);
		if(sz(A)>=sz(B)) rep(i,0,sz(A))rep(j,0,sz(B)) C[i+j]+=A[i]*B[j];
		else             rep(i,0,sz(B))rep(j,0,sz(A)) C[i+j]+=B[i]*A[j];
		return C;
	}
	static void mult_fft(FormalPowerSeries<tfps> &A, FormalPowerSeries<tfps> B){
		// A'=A*B in O(nlgn)
		const int nm=A.size()+B.size();
		const int lgi=ilog2(nm-1)+1;
		A.F.resize(1<<lgi,0);B.F.resize(1<<lgi,0);
		internal::fft(A.F,0);internal::fft(B.F,0);
		rep(i,0,sz(A))A[i]*=B[i];
		internal::fft(A.F,1);
		A.trunc(nm);
	}
	
	static void mult(FormalPowerSeries<tfps> &A, const FormalPowerSeries<tfps> &B){
		static const int Limit=40;
		if(min(sz(A),sz(B))<=Limit)A=mult_naive(A,B);
		else mult_fft(A,B);
	}
	
	static void scale(FormalPowerSeries<tfps> &F,tfps vl){
		// F*vl
		for(tfps &ac:F)ac*=vl;
	}
	static void add(FormalPowerSeries<tfps> &A, const FormalPowerSeries<tfps> &B, int sgn=1 ){
		// A + B*sgn ; I assume sgn\in\{-1,1\}
		A.F.resize(max(sz(A),sz(B)));
		if(sgn==1)      rep(i,0,min(sz(A),sz(B)))A[i]+=B[i];
		else if(sgn==-1)rep(i,0,min(sz(A),sz(B)))A[i]-=B[i];
		else assert(0);
	}
	
	
	static FormalPowerSeries<tfps> shift(const FormalPowerSeries<tfps> &F,int xi){
		// G=x^{xi}F
		FormalPowerSeries<tfps> G(max(1,sz(F)+xi),0);
		rep(i,0,sz(F))if(i+xi>=0&&i+xi<sz(G)) G[i+xi]=F[i];
		// idebug(F.F);idebug(G.F);
		return G;
	}
	
	// constructor
	constexpr FormalPowerSeries()noexcept:F({0}){};
	constexpr FormalPowerSeries(const std::vector<tfps> &Fp):F(Fp){};
	FormalPowerSeries(std::initializer_list<tfps> Fp):F(Fp){};
	constexpr FormalPowerSeries(int n,tfps vl=tfps(0)){F.resize(n,vl);}
	
	// iterator stuff
	using iterator=typename std::vector<mint>::iterator;
	using const_iterator=typename std::vector<mint>::const_iterator;
	iterator begin() { return F.begin(); } iterator end() { return F.end(); }
	const_iterator begin() const { return F.begin(); } const_iterator end() const { return F.end(); }
	
	// utilities
	FormalPowerSeries<tfps>& trunc(int n,bool elim_0=1){
		F.resize(max(1,min(sz(F),n)));
		if(elim_0)while(sz(F)>1&&F.back()==0)F.pop_back();
		return *this;
	}
	int size()const{return F.size();}
	tfps &operator[](int ind){F.resize(max(size(),ind+1),0);return F[ind];}
	tfps operator[](int ind)const {return (size()>ind)?F[ind]:tfps(0);}
	
	// basic operators
	FormalPowerSeries<tfps> &operator *=(const FormalPowerSeries<tfps>&B){mult(*this,B); return *this;}
	FormalPowerSeries<tfps> operator * (const FormalPowerSeries<tfps>&B)const{FormalPowerSeries<tfps> tmp=*this; return (tmp*=B);}
	FormalPowerSeries<tfps> &operator +=(const FormalPowerSeries<tfps>&B){add(*this,B); return *this;}
	FormalPowerSeries<tfps> operator + (const FormalPowerSeries<tfps>&B)const{FormalPowerSeries<tfps> tmp=*this; return (tmp+=B);}
	FormalPowerSeries<tfps> &operator -=(const FormalPowerSeries<tfps>&B){add(*this,B,-1); return *this;}
	FormalPowerSeries<tfps> operator - (const FormalPowerSeries<tfps>&B)const{FormalPowerSeries<tfps> tmp=*this; return (tmp-=B);}
	
	FormalPowerSeries<tfps> operator -()const{return FormalPowerSeries<tfps>({0}) - (*this);}
	
	//scalar operators
	FormalPowerSeries<tfps> &operator*=(const tfps &c){scale(*this,c);return *this;}
	FormalPowerSeries<tfps> operator * (const tfps &c)const{FormalPowerSeries<tfps> tmp=*this;return (tmp*=c);}
	
	// comparators
	bool operator ==(const FormalPowerSeries<tfps>&B)const{
		bool res=1;rep(i,0,max(size(),sz(B)))res=res&&(*this)[i]==B[i];
		return res;
	}
	bool operator !=(const FormalPowerSeries<tfps>&B)const{return !((*this)==B);}
	
	// shifting operators (x^{shf}*F-> F>>shf or F<<shf depending on shfs sign)
	FormalPowerSeries<tfps> operator<<(int shf)const{return shift(*this,shf);}
	FormalPowerSeries<tfps> &operator<<=(int shf){*this=(*this)<<shf;return *this;}
	FormalPowerSeries<tfps> operator>>(int shf)const{return shift(*this,-shf);}
	FormalPowerSeries<tfps> &operator>>=(int shf){*this=(*this)>>shf;return *this;}
	
	// square - derivative - integration - inverse - log - exp - sqrt
	
	FormalPowerSeries<tfps>& square(){
		const int Limit=70;
		if(size()<=Limit)(*this)=mult_naive(*this,*this);
		else {
			const int nm=size()*2, lgi=ilog2(nm-1)+1;
			F.resize(1<<lgi,0);
			internal::fft(F,0);
			for(tfps&ac:F)ac*=ac;
			internal::fft(F,1);
			trunc(nm);
		}
		return *this;
	}
	
	FormalPowerSeries<tfps> deriv()const{
		FormalPowerSeries<tfps> G(max(1,sz(F)-1),0);// G=D(F)
		rep(i,1,sz(F))G[i-1]=F[i]*tfps(i);
		return G;
	}
	FormalPowerSeries<tfps> integ()const{
		FormalPowerSeries<tfps> G(sz(F)+1);//D(G)=F
		rep(i,0,sz(F))G[i+1]=F[i]/tfps(i+1);
		return G;
	}
	
	FormalPowerSeries<tfps> inv(int n) const {
		assert(F[0]);// G*F=1
		FormalPowerSeries<tfps> G={F[0].inv()};//mpow(F[0],mod-2,mod)
		for(int e=2;e<2*n;e<<=1){
			FormalPowerSeries<tfps> ac=G;
			ac*= FormalPowerSeries({F.begin(),F.begin()+min(sz(F),e)}); // gives a ~/2 speedup
			for(tfps &act:ac)act=-act;
			ac+={2}; G*=ac;
			G.trunc(e);
		}
		return G.trunc(n);
	}
	FormalPowerSeries<tfps> log(int n) const {
		//first n terms of G=ln(F)=integ(D(F)/F)
		assert(F[0]==1);
		return (inv(n)*deriv()).trunc(n-1).integ();
	}

	FormalPowerSeries<tfps> exp(int n) const {
		assert(!F[0]);// first n terms of G=exp(F) 
		FormalPowerSeries<tfps> G={1},ac;
		for(int e=2;e<n*2;e<<=1){
			ac=G.log(e);ac.F.resize(max(sz(ac),min(e,sz(F))),0);
			
			rep(i,0,sz(ac))ac[i]=(i<sz(F)?F[i]:mint(0))- ac[i];
			ac+={1}; (G*=ac).trunc(e);
		}
		return G.trunc(n);
	}
	FormalPowerSeries<tfps> pow(lint e,int n)const{
		// first n terms of G=F^e O(nlgn)
		assert(e>=0); if(!e)return {1};
		int xi=0;while(xi<sz(F)&&F[xi]==0)++xi;
		if( xi>=sz(F) || lint(xi) > lint(n)/e )return {0};
		
		// alp could go wrong if its not a integer modulo smthng
		tfps alp=F[xi].pow(e%(mod-1)),ainv=F[xi].inv();
		
		int rx=xi?xi*int(e):0; // F^i=exp(i*log(F)) ; shifts and alp,ainv  for making F[0]=1
		return ( ( ( ((*this)>>xi)*ainv).trunc(n-xi).log(n-rx)*tfps(e) 
				).exp(n-rx) >> rx )*alp;
	}
	FormalPowerSeries<tfps> sqrt(int n)const {
		// this impl only works on numbers modulo a prime
		
		// first n terms of G^2=F O(nlgn)
		const tfps i2 = tfps(2).inv();//mpow(2,mod-2,mod);
		
		int xi=0;while(xi<sz(F)&&F[xi]==0)++xi;
		if( xi>=sz(F) )return {0};
		
		if(xi&1)return {};
		int vl=Tonelli_Shanks(F[xi]);if(vl==-1)return {};
		FormalPowerSeries<tfps> G={vl},ac;//sqrt{F[0]}
		//I assume Tonelli_Shanks returns -1 when it's impossible
		
		FormalPowerSeries<tfps> H=(*this)>>xi;H.trunc(n-xi);
		
		for(int e=2;e<(n-xi)*2;e<<=1){
			ac=G.inv(e)*H;ac.trunc(e);
			G+=ac; G*=i2; G.trunc(e);
		}
		return (G<<(xi/2)).trunc(n);
	}
	
};typedef FormalPowerSeries<mint> fps;

auto take_time=[&](){return std::chrono::high_resolution_clock::now();};
auto get_durat=[&](auto start){ return std::chrono::duration_cast<std::chrono::nanoseconds>(take_time() - start).count(); };
//std::chrono::milliseconds
//std::chrono::microseconds
//std::chrono::nanoseconds
void multiplication() {
	const int m=1e4;
	// const int m=100;
	const int u=30,d=120;
	rep(i,u,d){
		fps A(i);rep(j,0,i)A[j]=rng_64();
		fps B(m);rep(j,0,m)B[j]=rng_64();
		
		fps Ap=A,Bp=B;
		
		auto st1=take_time();
		fps::mult_fft(A,B);
		// A.square(0);
		auto t1=get_durat(st1);
		auto st2=take_time();
		// Ap.square(1);
		fps::mult_naive(Ap,Bp);
		auto t2=get_durat(st2);
		debug(t1<t2,i,t1,t2);
	}
}
void expon(){
	const int n=1e1,m_pow=10;
	fps F(n);rep(i,0,n)F[i]=rng_64();
	rep(e,0,m_pow){
		// int ex=(1<<e)-1;
		int ex=2;
		auto st1=take_time();
		fps G1=F.pow(ex,n);
		auto t1=get_durat(st1);
		
		// debug("xd");
		
		auto st2=take_time();
		fps G2=F.pow2(ex,n);
		auto t2=get_durat(st2);
		
		idebug(G1);idebug(G2);
		assert(G1==G2);
		debug(t1,t2,e);
	}
	
}

int32_t main(){
	ios_base::sync_with_stdio(false);
    cin.tie(NULL);
	cout << setprecision(12) << fixed;

    // multiplication();
	expon();
	return 0;
}

