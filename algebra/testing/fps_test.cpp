// autor:Oscar Vargas Pabon
// control impl ends at line 500

#include <bits/stdc++.h>
#include<cassert>

typedef long long lint;

using namespace std;

#define debug(args...) { string _s = #args; replace(_s.begin(), _s.end(), ',', ' '); stringstream _ss(_s); istream_iterator<string> _it(_ss); raw_debug(_it, args);}
void raw_debug(istream_iterator<string> it) {cerr<<endl;assert(it==it);}
template<typename T, typename... Args>
void raw_debug(istream_iterator<string> it, T a, Args... args) { cerr <<"<"<< *it << "->" << a << "> "; raw_debug(++it, args...); }
#define idebug(v) {cout<<'['<<#v<<']';for(const auto &el:v)cout << ' ' << el; cout << endl;}
#define adebug(ar,n) {cout<<'['<<#ar<<']';for(int i=0;i<n;++i)cout << ' ' << ar[i]; cout << endl;}

#define rep(i,strt,end) for(int i = strt ; i !=int(end) ; (int(strt)<int(end))?++i:--i )
#define rall(vec) vec.rbegin(), vec.rend()
#define all(vec) vec.begin(), vec.end()
#define sz(vec) int(vec.size())
#define pb push_back
#define pob pop_back
#define pf push_front
#define pof pop_front

const auto rng_seed=chrono::steady_clock::now().time_since_epoch().count();
mt19937_64 rng_64( rng_seed );
constexpr int ilog2( int num ) { return 8*sizeof(int) - __builtin_clz( num ) - 1; }
constexpr int mpow(int x,int e,int m){int res=1;while(e){if(e&1)res=(res*1ll*x)%m;e>>=1;x=(x*1ll*x)%m;}return res;}

namespace control{

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

I assume from my template:
everything     :: #define rep(i,strt,end) for(int i = strt ; i !=int(end) ; (int(strt)<int(end))?++i:--i )
everything     :: #define sz(vec) int(vec.size())
p_trunc        :: #define pob pop_back
p_mult,p_square:: int ilog2( int num ) { return 8*sizeof(int) - __builtin_clz( num ) - 1; }
modInverse     :: int mpow(int x,int e,int m){int res=1;while(e){if(e&1)res=(res*1ll*x)%m;e>>=1;x=(x*1ll*x)%m;}return res;}	
TonelliShanks  :: mt19937_64 rng_64( chrono::steady_clock::now().time_since_epoch().count() );
*/
struct mint{
	const static int mod= 998244353; static_assert(mod>0);
	
	int vl;
	constexpr mint()noexcept:vl(0){};
	constexpr mint( int v)noexcept:vl(v>=0?v%mod:(v%mod)+mod){};
	constexpr mint(lint v)noexcept:vl(v>=0?v%mod:(v%mod)+mod){};
	constexpr mint(unsigned long long v)noexcept:vl(v%mod){};
	
	mint &operator +=(const mint &ot){ vl+=ot.vl; if(vl>=mod)vl-=mod; return *this; }
	mint  operator + (const mint &ot)const{ return mint(*this)+=ot; }
	mint &operator -=(const mint &ot){ vl-=ot.vl; if(vl<0)vl+=mod; return *this; }
	mint  operator - (const mint &ot)const{ return mint(*this)-=ot; }
	mint &operator *=(const mint &ot){ vl=(vl*1ll*ot.vl)%mod; return *this; }
	mint  operator * (const mint &ot)const{ return mint(*this)*=ot; }
	mint &operator /=(const mint &ot){ (*this)*=ot.inverse(); return *this; }
	mint  operator / (const mint &ot)const{ return mint(*this)/=ot; }
	
	mint inverse()const{return mint(mpow(vl,mod-2,mod));}//Fermats little theorem
	mint operator -()const {return mint(-vl);}
	mint pow(lint e)const{return mint(mpow(vl,e%(mod-1),mod));}
	
	bool operator ==(const mint &ot)const{return vl==ot.vl;} 
	bool operator ==(const  int &ot)const{return vl==ot;   }
	bool operator !=(const mint &ot)const{return vl!=ot.vl;}
	
	operator bool() const { return vl; }
	operator  int() const { return vl; }
	
	friend ostream &operator<<(ostream &os,const mint &ac){return os << ac.vl;}
	friend istream &operator>>(istream&is,mint &ac){int v;cin>>v;ac=mint(v);return is;}	
}; const int mod=mint::mod;

typedef mint tfps; typedef vector<tfps> fps; // fps renaming types

/* START OF NTT */
namespace internal { // taken from atcoder

int countr_zero(unsigned int n) { return __builtin_ctz(n); }
constexpr int countr_zero_constexpr(unsigned int n) { int x = 0; while (!(n & (1 << x))) x++; return x; }

struct fft_info {
	static const int g =3;//mpow(3,119,mod); // primitive root 
	
    static constexpr int rank2 = countr_zero_constexpr(mint::mod - 1);
    std::array<mint, rank2 + 1> root;   // root[i]^(2^i) == 1
    std::array<mint, rank2 + 1> iroot;  // root[i] * iroot[i] == 1

    std::array<mint, std::max(0, rank2 - 2 + 1)> rate2;
    std::array<mint, std::max(0, rank2 - 2 + 1)> irate2;

    std::array<mint, std::max(0, rank2 - 3 + 1)> rate3;
    std::array<mint, std::max(0, rank2 - 3 + 1)> irate3;

    fft_info() {
        root[rank2] = mint(g).pow((mint::mod - 1) >> rank2);
        iroot[rank2] = root[rank2].inverse();
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


void butterfly(fps& a) {
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
                    auto mod2 = 1ULL * mint::mod * mint::mod;
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
                        (unsigned long long)((unsigned int)(l.vl - r.vl) + mint::mod) *
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
                        mint((mint::mod + a2 - a3) * iimag.vl).vl;

                    a[i + offset] = a0 + a1 + a2 + a3;
                    a[i + offset + 1 * p] =
                        (a0 + (mint::mod - a1) + a2na3iimag) * irot.vl;
                    a[i + offset + 2 * p] =
                        (a0 + a1 + (mint::mod - a2) + (mint::mod - a3)) *
                        irot2.vl;
                    a[i + offset + 3 * p] =
                        (a0 + (mint::mod - a1) + (mint::mod - a2na3iimag)) *
                        irot3.vl;
                }
                if (s + 1 != (1 << (len - 2)))
                    irot *= info.irate3[countr_zero(~(unsigned int)(s))];
            }
            len -= 2;
        }
    }
}

} //end internal namespace
void fft(fps&a,bool invert){
	if(invert){
		internal::butterfly_inv(a);
        mint n_1 = mint(sz(a)).inverse();
        for (mint & x : a)x*=n_1;
	} else internal::butterfly(a);
}
/* End of NTT */

void p_trunc(fps &F, int n,bool elim_0=1){
	F.resize(max(1,min(sz(F),n)));
	if(elim_0)while(sz(F)>1&&F.back()==0)F.pob();
}
fps p_mult_naive(const fps&A,const fps &B){
	fps C(sz(A)+sz(B),0); //C[i+j]=(C[i+j]+A[i]*1ll*B[j])%mod;
	if(sz(A)>=sz(B)) rep(i,0,sz(A))rep(j,0,sz(B)) C[i+j]+=A[i]*B[j];
	else             rep(i,0,sz(B))rep(j,0,sz(A)) C[i+j]+=B[i]*A[j];
	return C;
}
/*fps p_mult_sparse(const fps&A,const fps &B){
	fps C(sz(A)+sz(B),0);//deprecate this shit, its useless
	vector<int> tmp;rep(i,0,sz(A))if(A[i])tmp.push_back(i);
	rep(i,0,sz(B))if(B[i]){
		// for(int j:tmp)C[i+j]=(C[i+j]+ B[i]*1ll*A[j])%mod;
		for(int j:tmp)C[i+j]+=B[i]*A[j];
	}
	return C;
}*/
void p_mult_fft(fps &A, fps B){
	// A'=A*B in O(nlgn)
	int nm=A.size()+B.size();
	int lgi=ilog2(nm-1)+1;A.resize(1<<lgi,0);B.resize(1<<lgi,0);
	fft(A,0);fft(B,0);
	// rep(i,0,sz(A))A[i]=(A[i]*1ll*B[i])%mod;
	rep(i,0,sz(A))A[i]*=B[i];
	fft(A,1);
	p_trunc(A,nm);
}
void p_mult(fps &A, const fps &B){
	/*auto cnt_zeros=[](const fps &F){int cnt=0;for(tfps ac:F){cnt+=ac==0;if(cnt>=Limit)return false;}return true;};
	else if(cnt_zeros(A)||cnt_zeros(B)) A=p_mult_sparse(A,B);*/
	static const int Limit=50;
	if(min(sz(A),sz(B))<=Limit)A=p_mult_naive(A,B);
	else p_mult_fft(A,B);
}
void p_square_fft(fps &A){
	// A'=A*A in O(nlgn)
	int nm=A.size()*2;
	int lgi=ilog2(nm-1)+1;A.resize(1<<lgi);
	fft(A,0);
	// for(int &ac:A)ac=(ac*1ll*ac)%mod;
	for(tfps &ac:A)ac*=ac;
	fft(A,1);
	p_trunc(A,nm);
}
void p_square(fps &A){
	static const int Limit=50;
	if(sz(A)<=Limit)A=p_mult_naive(A,A);
	else p_square_fft(A);
}

void p_add(fps &F, tfps vl, int xi=0){
	// F+ vl*x^{xi}
	if(sz(F)<=xi)F.resize(xi+1,0);
	F[xi]+=vl;
}
void p_scale(fps &F,tfps vl){
	// F*vl
	// rep(i,0,sz(F))F[i]=(F[i]*1ll*vl)%mod;
	for(tfps &ac:F)ac*=vl;
}
void p_add(fps &A, const fps &B, int sgn=1 ){
	// A + B*sgn ; I assume sgn\in\{-1,1\}
	A.resize(max(sz(A),sz(B)));
	// rep(i,0,min(sz(A),sz(B)))A[i]=(A[i]+B[i]*sgn)%mod;
	// rep(i,0,min(sz(A),sz(B)))if(A[i]<0)A[i]+=mod;
	if(sgn==1)      rep(i,0,min(sz(A),sz(B)))A[i]+=B[i];
	else if(sgn==-1)rep(i,0,min(sz(A),sz(B)))A[i]-=B[i];
	else assert(0);
}
fps p_deriv(const fps &F){
	fps G(max(1,sz(F)-1));// G=D(F)
	// rep(i,1,sz(F))G[i-1]=(F[i]*1ll*i)%mod;
	rep(i,1,sz(F))G[i-1]=F[i]*tfps(i);
	return G;
}
fps p_inte(const fps &F){
	fps G(sz(F)+1);//D(G)=F
	// rep(i,0,sz(F))G[i+1]=(F[i]*1ll*mpow(i+1,mod-2,mod))%mod;
	rep(i,0,sz(F))G[i+1]=F[i]/tfps(i+1);
	return G;
}
fps p_shift(const fps &F,int xi){
	// G=x^{xi}F
	fps G(max(1,sz(F)+xi),0);
	rep(i,0,sz(F))if(i+xi>=0&&i+xi<sz(G)) G[i+xi]=F[i];
	return G;
}

fps p_inv(const fps &F,int n){
	assert(F[0]);// G*F=1
	fps G={F[0].inverse()};//mpow(F[0],mod-2,mod)
	for(int e=2;e<2*n;e<<=1){
		fps ac=G;// gives a ~/2 speedup
		p_mult(ac,{F.begin(),F.begin()+min(sz(F),e)});
		// rep(i,0,sz(ac))ac[i]=(mod-ac[i])%mod;
		for(tfps &act:ac)act=-act;
		
		
		p_add(ac,2); p_mult(G,ac);
		p_trunc(G,e,0);
	}
	p_trunc(G,n);
	return G;
}
fps p_log(const fps &F,int n){
	//first n terms of G=ln(F)=inte(D(F)/F)
	assert(F[0]==1);
	fps G=p_inv(F,n),dF=p_deriv(F);
	p_mult(G,dF); p_trunc(G,n-1);
	return p_inte(G);
}

fps p_exp(const fps &F,int n){
	assert(!F[0]);// first n terms of G=exp(F) 
	fps G={1},ac;
	for(int e=2;e<n*2;e<<=1){
		ac=p_log(G,e); ac.resize(max(sz(ac),min(e,sz(F))),0);
		rep(i,0,sz(ac))ac[i]=(i<sz(F)?F[i]:mint(0))- ac[i];
		// rep(i,0,sz(ac))ac[i]=(
			// (i<sz(F)?F[i]:0)- ac[i])%mod;
		// for some reason working with negatives
		// destroys something (NTT???)
		// rep(i,0,sz(ac))if(ac[i]<0)ac[i]+=mod;
		p_add(ac,1);
		
		p_mult(G,ac); p_trunc(G,e,0);
	}
	p_trunc(G,n);
	return G;
}
fps p_lgipow(const fps &F,lint e,int n){
	// first n terms of G=F^e O(nlgn)
	if(!e)return {1};
	int xi=0;while(xi<sz(F)&&F[xi]==0)++xi;
	if( xi>=sz(F) || xi > n/e )return {0};
	
	tfps alp=mpow(F[xi],e%(mod-1),mod),ainv=mpow(F[xi],mod-2,mod);
	fps H(F.begin()+xi,F.end());p_scale(H,ainv);
	p_trunc(H,n-xi);
	
	int rx=xi*int(e);
	H=p_log(H,n-rx);
	p_scale(H,tfps(e));
	H=p_exp(H,n-rx);
	
	fps G(n,0); H.resize(n-rx);
	// rep(i,rx,n)G[i]=(H[i-rx]*1ll*alp)%mod;
	rep(i,rx,n)G[i]=H[i-rx]*alp;
	return G;
}
fps p_binpow( fps F, lint e, int n){
	// first n terms of G=F^e O(nlgnlge)
	fps G={1};while(e){
		if(e&1ll)p_mult(G,F);
		e>>=1;p_square(F);
		p_trunc(G,n);p_trunc(F,n);
	}
	return G;
}
fps p_pow(const fps&F,lint e,int n){
	if(e<(1ll<<2))return p_binpow(F,e,n);
	else return p_lgipow(F,e,n);
}


int Tonelli_Shanks(mint a) {
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
		// int error = (((mpow(a, mod - 2, mod) *1ll* x) % mod) *1ll* x) % mod;
		mint error= mint(mpow(a,mod-2,mod))*x*x;
		if (mpow(error, 1 << (Q - shift), mod) != 1) x *= b;
		b *=b;
		++shift;
	}
	return x;
}


fps p_sqrt( const fps &F,int n){
	// first n terms of G^2=F O(nlgn)
	const mint i2 = mint(2).inverse();//mpow(2,mod-2,mod);
	
	int xi=0;while(xi<sz(F)&&F[xi]==0)++xi;
	if( xi>=sz(F) )return {0};
	
	if(xi&1)return {};
	int vl=Tonelli_Shanks(F[xi]);if(vl==-1)return {};
	fps G={vl},ac;//sqrt{F[0]}
	//I assume Tonelli_Shanks returns -1 when it's impossible
	
	fps H(F.begin()+xi,F.end()); p_trunc(H,n-xi);
	
	for(int e=2;e<(n-xi)*2;e<<=1){
		ac=p_inv(G,e);
		p_mult(ac,H);p_trunc(ac,e,0);
		p_add(G,ac);p_scale(G,i2);
		p_trunc(G,e);
	}
	swap(H,G);
	G.resize(sz(G)+xi/2);
	rep(i,0,xi/2)G[i]=0;
	rep(i,xi/2,sz(G))G[i]=H[i-xi/2];
	return G;
}


}// end namespace control
#undef sz
/******* end of control impl *******/
template<typename tpow> constexpr tpow mpow(tpow x,lint e,tpow m){tpow res=1;while(e){if(e&1ll)res=(res*1ll*x)%m;e>>=1;x=(x*1ll*x)%m;}return res;}

const int mod = 998244353;/*
Author: Oscar Vargas Pabon

It is though to work modulo primes, so inv (.inv,/,/=) and pow may not
	work properly otherwise.

I assume from my template :
inv :: int mpow(int x,int e,int m){int res=1;while(e){if(e&1)res=(res*1ll*x)%m;e>>=1;x=(x*1ll*x)%m;}return res;}	

Tested in testing/test_alghelp.cpp and in fft/ntt stuff
*/

template<lint raw_m,typename tint=int,bool arbi_ntt=0>
struct modulo_int{
	constexpr static tint m=raw_m; static_assert(m>0);
	constexpr static tint mod(){return m;}
	constexpr static bool arbitrary_ntt(){return arbi_ntt;}
	
	tint vl;
	inline constexpr modulo_int()noexcept:vl(0){};
	inline constexpr modulo_int( int v)noexcept:vl(v>=0?(v<m?v:v%m):(v+m>=0?v+m:(v%m)+m)){};
	inline constexpr modulo_int(lint v)noexcept:vl(v>=0?(v<m?v:v%m):(v+m>=0?v+m:(v%m)+m)){};
	inline constexpr modulo_int(unsigned long long v)noexcept:vl(v<m?v:v%m){};
	
	inline constexpr modulo_int &operator +=(const modulo_int &ot){ vl+=ot.vl; if(vl>=m)vl-=m; return *this; }
	inline constexpr modulo_int  operator + (const modulo_int &ot)const{ return modulo_int(*this)+=ot; }
	inline constexpr modulo_int &operator -=(const modulo_int &ot){ vl=(vl>=ot.vl)?vl-ot.vl:m-ot.vl+vl; return *this; }
	inline constexpr modulo_int  operator - (const modulo_int &ot)const{ return modulo_int(*this)-=ot; }
	inline constexpr modulo_int &operator *=(const modulo_int &ot){ vl=(vl*1ll*ot.vl)%m; return *this; }
	inline constexpr modulo_int  operator * (const modulo_int &ot)const{ return modulo_int(*this)*=ot; }
	inline constexpr modulo_int &operator /=(const modulo_int &ot){ (*this)*=ot.inv(); return *this; }
	inline constexpr modulo_int  operator / (const modulo_int &ot)const{ return modulo_int(*this)/=ot; }
	
	inline constexpr modulo_int inv()const{return modulo_int(vl).pow(m-2);}//Fermats little theorem
	inline constexpr modulo_int operator -()const {return modulo_int(-vl);}
	inline constexpr modulo_int pow(lint e)const{return modulo_int(mpow<tint>(vl,e%(m-1),m));}
	
	inline constexpr bool operator ==(const modulo_int &ot)const{return vl==ot.vl;} 
	inline constexpr bool operator ==(const  int &ot)const{return vl==ot;   }
	inline constexpr bool operator !=(const modulo_int &ot)const{return vl!=ot.vl;}
	
	inline constexpr operator bool() const { return vl; }
	inline constexpr operator  int() const { return vl; }
	inline constexpr operator long long() const { return vl; }
	
	friend ostream &operator<<(ostream &os,const modulo_int &ac){return os << ac.vl;}
	friend istream &operator>>(istream&is,modulo_int &ac){int v;is>>v;ac=modulo_int(v);return is;}	
}; typedef modulo_int<mod> mint;

/*
Author: Oscar Vargas Pabon
NTT made by yours truly under the help of several books and
	the almighty idea behind the recent progress on 
	https://codeforces.com/blog/entry/151162
previously, NTT was taken from the atcoder library
	https://github.com/atcoder/ac-library/
Operations on FPST based on code by MarcosK, other people and multiple blogs all around the place
	https://codeforces.com/contest/438/submission/340901913
	https://cp-algorithms.com/algebra/polynomial.html#inverse-series_1
	https://codeforces.com/blog/entry/56422
	https://codeforces.com/blog/entry/12513 - problem E
Tested in
	https://judge.yosupo.jp/problem/inv_of_formal_power_series
	https://judge.yosupo.jp/problem/exp_of_formal_power_series
	https://judge.yosupo.jp/problem/log_of_formal_power_series
	https://judge.yosupo.jp/problem/pow_of_formal_power_series
	https://judge.yosupo.jp/problem/sqrt_of_formal_power_series
Though this last version was tested in /testing/fps_test.cpp
 
Arbitrary modulus multiplication tested/developed in 
	https://atcoder.jp/contests/arc215/submissions/73551697
	with help of https://nyaannyaan.github.io/library/ntt/arbitrary-ntt.hpp
 
I assume mint from algebra/modulo_int.cpp
Note that in this specific version, the ntt somehow gets slower with montgomery_space (its
	safe but slower smhw).
 
Remember 'primitive_root(mod)' will return a primitive root of the modulo in time
	$O(\sqrt{mod}+R log^2 mod)$ usefull to get the constant
 
Dependencies I assume
everything     :: #define rep(i,strt,end) for(int i = strt ; i !=int(end) ; (int(strt)<int(end))?++i:--i )
mult,square    :: int ilog2( int num ) { return 8*sizeof(int) - __builtin_clz( num ) - 1; }
TonelliShanks  :: mt19937_64 rng_64( chrono::steady_clock::now().time_since_epoch().count() );
sqrt           :: TonelliShanks, tfps(?).inv()
pow            :: tfps(?).inv(), tfps(?).pow()
 
Notes: Some undefined behaviour may ocur whenever F[ind1]=F[ind2] where F.size()<=ind2 as
		F[ind1] can be evaluated first and the resize in F[ind2] may 'overwrite' the
		reference previously established
		
Notes: (1<<59)*27+1 is a 64-bit prime (IntroMathComputational-Shoup pg 484)
*/ namespace internal {
constexpr int countr_zero_constexpr(unsigned int n) { int x = 0; while (!(n & (1 << x))) x++; return x; }
constexpr int primitive_root(int m){
	//tests for g such that $\forall_{p|(md-1)} g^{(md-1)/p}=1(mod md)$
	int dec[32]={},dind=0, md=m-1;
	for(int i=2;i*i<=md;++i)if(md%i==0){
		dec[dind++]=i; while(md%i==0)md/=i;
	} if(md>1)dec[dind++]=md;
	
	bool fnd=0; int pr=1;while(!fnd){
		++pr; fnd=1; // by properties, this will always end
		for(int i=0;i<dind&&fnd;++i)fnd=mpow<int>(pr,(m-1)/dec[i],m)!=1;
	} return pr;// std::cerr <<pr << " _ primitive root" << endl;
} template<typename tfps,int rank> constexpr tfps fft_root(){
    constexpr tfps c_root=primitive_root(tfps::mod());
    constexpr int rank2=countr_zero_constexpr(tfps::mod()-1);
    if constexpr(rank>=rank2)
        return c_root.pow((tfps::mod()-1)>>rank);
    else return fft_root<tfps,rank+1>()*fft_root<tfps,rank+1>();
} template<typename tfps,int rank> constexpr tfps fft_iroot(){
    constexpr int rank2=countr_zero_constexpr(tfps::mod()-1);
    if constexpr(rank>=rank2)
        return fft_root<tfps,rank>().inv();
    else return fft_iroot<tfps,rank+1>()*fft_iroot<tfps,rank+1>();
} template<typename tfps,int n> void butterfly_rec(tfps*a){
    if constexpr(n<=1)return;
    constexpr int e=ilog2(n),m=n/2;
    constexpr tfps wlen=fft_root<tfps,e>();
    tfps w=1; for(int i=0;i<m;++i){
        const tfps u=a[i],v=a[i+m];
        a[i  ]=u+v; a[i+m]=(u-v)*w;
        w*=wlen;
    }butterfly_rec<tfps,m>(a);butterfly_rec<tfps,m>(a+m);
}template<typename tfps,int n=(1<<countr_zero_constexpr(tfps::mod()-1))>
void butterfly(std::vector<tfps>&a){
    if constexpr(n<=1)return;
    if (n==int(a.size()))butterfly_rec<tfps,n>(a.data());
    else butterfly<tfps,n/2>(a);
} template<typename tfps,int n> void ibutterfly_rec(tfps*a){
    if constexpr(n<=1)return;
    constexpr int e=ilog2(n),m=n/2;
    constexpr tfps wlen=fft_iroot<tfps,e>();
    ibutterfly_rec<tfps,m>(a); ibutterfly_rec<tfps,m>(a+m);
    tfps w=1;for(int i=0;i<m;++i){
        const tfps u=a[i],v=w*a[i+m];
        a[i  ]=u+v; a[i+m]=u-v;
        w*=wlen; }
}template<typename tfps,int n=(1<<countr_zero_constexpr(tfps::mod()-1))>
void ibutterfly(std::vector<tfps>&a){
    if constexpr(n<=1)return;
    if (n==int(a.size()))ibutterfly_rec<tfps,n>(a.data());
    else ibutterfly<tfps,n/2>(a);
} template<typename tfps> void fft(vector<tfps>&a,bool invert){
    if(invert){ internal::ibutterfly(a);
        tfps n_1 = tfps(int(a.size())).inv();
        for (tfps&x:a)x*=n_1;
	}else internal::butterfly(a);
} template<typename tfps> void transposed_fft(std::vector<tfps> &A,bool invert){
	if (!invert) internal::ibutterfly<tfps>(A);
	reverse(A.begin() + 1, A.end()); if(invert){
		internal::butterfly<tfps>(A);
		tfps n_1=tfps(A.size()).inv();
		for (tfps &x: A) x *= n_1; }
} template<typename tfps,const bool square=0>
void mult_arbitrary(std::vector<tfps>&A,const std::vector<tfps> &B={} ){
	// Im guiding miself on nyaan's library for this
	// https://nyaannyaan.github.io/library/ntt/arbitrary-ntt.hpp
	using u128 = __uint128_t;
	const int nm=A.size()+B.size(),lgi=ilog2(nm-1)+1,n=1<<lgi,bsz=B.size();
	A.resize(n,0); // IntroMathComputational-Shoup pg 484 __ they get overflow
	// constexpr int m0 = (1<<30)*3+1, m1 = (1<<28)*13+1, m2 = (1<<27)*29+1;
	constexpr int m0 = 167772161, m1 = 469762049, m2 = 754974721;
	
	std::vector<modulo_int<m0>> A0(n);
	for(int i=0;i<n;++i)A0[i]=lint(A[i]);
	internal::fft<modulo_int<m0>>(A0,0);
	
	if constexpr(square){
		rep(i,0,1<<lgi)A0[i]*=A0[i];
	} else{
		std::vector<modulo_int<m0>> B0(n);
		for(int i=0;i<n;++i)B0[i]=(i>=bsz)?0:lint(B[i]);
		internal::fft<modulo_int<m0>>(B0,0);
		rep(i,0,1<<lgi)A0[i]*=B0[i];	
	}internal::fft<modulo_int<m0>>(A0,1);
	
	if constexpr ( tfps::mod()<m0 && tfps::mod()*1ll*tfps::mod() < m0 ){
		for(int i=0;i<n;++i)A[i]=int(A0[i]);
		return;
	} std::vector<modulo_int<m1>> A1(1<<lgi);
	for(int i=0;i<n;++i)A1[i]=lint(A[i]);
	internal::fft<modulo_int<m1>>(A1,0);
	
	if constexpr(square){
		rep(i,0,1<<lgi)A1[i]*=A1[i];
	} else {
		std::vector<modulo_int<m1>> B1(1<<lgi);
		for(int i=0;i<n;++i)B1[i]=(i>=bsz)?0:lint(B[i]);
		internal::fft<modulo_int<m1>>(B1,0);
		rep(i,0,1<<lgi)A1[i]*=B1[i];
	}internal::fft<modulo_int<m1>>(A1,1);
		
		
	constexpr bool only_2=0,cond_2=tfps::mod()<max<int>(m0,m1) && tfps::mod()*1ll*tfps::mod() < m0*1ll*m1;
	if constexpr ( only_2 || cond_2 ){
		constexpr int i01=modulo_int<m0>(m1).inv(),i10=modulo_int<m1>(m0).inv();
		constexpr lint m01=m0*1ll*m1;
		for(int i=0;i<n;++i){
			A[i]= lint( (
			        int(A0[i])*u128( i01*1ll*m1 ) +
				    int(A1[i])*u128( i10*1ll*m0 ) ) %m01 );
		} return;
	} std::vector<modulo_int<m2>> A2(1<<lgi);
	for(int i=0;i<n;++i)A2[i]=lint(A[i]);
	internal::fft<modulo_int<m2>>(A2,0);
	if constexpr(square){
		rep(i,0,1<<lgi)A2[i]*=A2[i];
	} else {
		std::vector<modulo_int<m2>> B2(1<<lgi);
		for(int i=0;i<n;++i)B2[i]=(i>=bsz)?0:lint(B[i]);
		internal::fft<modulo_int<m2>>(B2,0);
		rep(i,0,1<<lgi)A2[i]*=B2[i];
	} internal::fft<modulo_int<m2>>(A2,1);
	
	constexpr int i12=modulo_int<m0>(m1*1ll*m2).inv(),
				  i02=modulo_int<m1>(m0*1ll*m2).inv(),
				  i01=modulo_int<m2>(m0*1ll*m1).inv();
	constexpr u128 m012=u128(m0)*m1*m2;
	for(int i=0;i<n;++i){
		A[i]=lint( ( (
			  int(A0[i])*u128( (m1*1ll*m2*u128(i12))%m012 ) +
			  int(A1[i])*u128( (m0*1ll*m2*u128(i02))%m012 ) +
			  int(A2[i])*u128( (m0*1ll*m1*u128(i01))%m012 ) )%m012 )%tfps::mod() );
	}
} template<typename tfps> int Tonelli_Shanks(tfps a) {
	// usado por sqrt cuando trabajo con enteros modulo algo
	//plagiado epicamente de https://judge.yosupo.jp/submission/270105
	const int mod=tfps::mod();
	if (int(a) < 2) return a;
	if (mpow<int>(a, (mod - 1) / 2, mod) != 1) return -1;
	if (mod % 4 == 3) return mpow<int>(a, (mod + 1) / 4, mod);
 
	tfps b = 3; if (mod != 998244353) {
		while (mpow<int>(b, (mod - 1) / 2, mod) == 1) {
			b=tfps(int(rng_64()%(mod-3)) + 2);
		}
	} int q = mod - 1,Q = 0;
	while ( !(q&1) ) Q++, q /= 2;
 
	tfps x = mpow<int>(a, (q + 1) / 2, mod);
	b = mpow<int>(b, q, mod);
 
	int shift = 2; while ( x*x != a) {
		tfps error= tfps(mpow<int>(a,mod-2,mod))*x*x;
		if (mpow<int>(error, 1 << (Q - shift), mod) != 1) x *= b;
		b *=b; ++shift;
	} return x; }
} //end internal namespace
template<typename tfps> struct FormalPowerSeries{
	std::vector<tfps> F;
	static FormalPowerSeries<tfps> mult_naive(const FormalPowerSeries<tfps>&A,const FormalPowerSeries<tfps>&B){
		FormalPowerSeries<tfps> C(A.size()+B.size(),0);
		if(A.size()>=B.size())rep(i,0,A.size())rep(j,0,B.size())C[i+j]+=A[i]*B[j];
		else                  rep(i,0,B.size())rep(j,0,A.size())C[i+j]+=B[i]*A[j];
		return C;
	} static void mult_fft(FormalPowerSeries<tfps> &A, FormalPowerSeries<tfps> B){
		// A'=A*B in O(nlgn)
		const int nm=A.size()+B.size();
		const int lgi=ilog2(nm-1)+1;
		A.F.resize(1<<lgi,0);B.F.resize(1<<lgi,0);
		internal::fft<tfps>(A.F,0);internal::fft<tfps>(B.F,0);
		rep(i,0,A.size())A[i]*=B[i];
		internal::fft<tfps>(A.F,1);
	} static void mult(FormalPowerSeries<tfps> &A, const FormalPowerSeries<tfps> &B){
		const int nm=A.F.size()+B.F.size();
		static const int Limit=20;
		if(min<int>(A.size(),B.size())<=Limit)A=mult_naive(A,B);
		else {
			if constexpr( tfps::arbitrary_ntt() )
				internal::mult_arbitrary<tfps>(A.F,B.F);
			else mult_fft(A,B);
		} A.trunc(nm);
	} static void scale(FormalPowerSeries<tfps> &F,tfps vl){
		for(tfps &ac:F)ac*=vl; // F*vl
	} static void add(FormalPowerSeries<tfps> &A, const FormalPowerSeries<tfps> &B, int sgn=1 ){
		A.F.resize(max<int>(A.size(),B.size())); // A + B*sgn ; I assume sgn\in\{-1,1\}
		if(sgn==1)      rep(i,0,min<int>(A.size(),B.size()))A[i]+=B[i];
		else if(sgn==-1)rep(i,0,min<int>(A.size(),B.size()))A[i]-=B[i];
		else assert(0);
	} static FormalPowerSeries<tfps> shift(const FormalPowerSeries<tfps> &F,int xi){
		FormalPowerSeries<tfps> G(max<int>(1,F.size()+xi),0); // G=x^{xi}F
		rep(i,0,F.size())if(i+xi>=0&&i+xi<G.size()) G[i+xi]=F[i];
		return G;
	} // ########################### constructor ##################################
	inline constexpr FormalPowerSeries()noexcept:F({0}){};
	inline constexpr FormalPowerSeries(const std::vector<tfps> &Fp):F(Fp){};
	inline FormalPowerSeries(std::initializer_list<tfps> Fp):F(Fp){};
	inline constexpr FormalPowerSeries(int n,tfps vl=tfps(0)){F.resize(n,vl);}
	// ############################ iterator stuff ################################
	using iterator=typename std::vector<tfps>::iterator;
	using const_iterator=typename std::vector<tfps>::const_iterator;
	iterator begin() { return F.begin(); } iterator end() { return F.end(); }
	const_iterator begin() const { return F.begin(); } const_iterator end() const { return F.end(); }
	// ############################## utilities ###################################
	FormalPowerSeries<tfps>& trunc(int n=-1,bool elim_0=1){
		if(n==-1)n=F.size();
		F.resize(max<int>(1,min<int>(F.size(),n)));
		if(elim_0)while(int(F.size())>1&&F.back()==0)F.pop_back();
		return *this;
	} int size()const{return F.size();}
	tfps &operator[](int ind){F.resize(max<int>(size(),ind+1),0);return F[ind];}
	tfps operator[](int ind)const {return (size()>ind)?F[ind]:tfps(0);}
	// ############################ basic operators ###############################
	inline FormalPowerSeries<tfps>&operator*=(const FormalPowerSeries<tfps>&B){mult(*this,B); return *this;}
	inline FormalPowerSeries<tfps> operator* (const FormalPowerSeries<tfps>&B)const{FormalPowerSeries<tfps> tmp=*this; return (tmp*=B);}
	inline FormalPowerSeries<tfps>&operator+=(const FormalPowerSeries<tfps>&B){add(*this,B); return *this;}
	inline FormalPowerSeries<tfps> operator+ (const FormalPowerSeries<tfps>&B)const{FormalPowerSeries<tfps> tmp=*this; return (tmp+=B);}
	inline FormalPowerSeries<tfps>&operator-=(const FormalPowerSeries<tfps>&B){add(*this,B,-1); return *this;}
	inline FormalPowerSeries<tfps> operator- (const FormalPowerSeries<tfps>&B)const{FormalPowerSeries<tfps> tmp=*this; return (tmp-=B);}
	inline FormalPowerSeries<tfps> operator% (const FormalPowerSeries<tfps>&B)const{return ((*this)-B*euc_div(B)).trunc(B.size()); }
	inline FormalPowerSeries<tfps>&operator%=(const FormalPowerSeries<tfps>&B)const{ (*this)=(*this)%B; return *this;}
	inline FormalPowerSeries<tfps> operator-()const{return FormalPowerSeries<tfps>({0}) - (*this);}
	//############################### scalar operators ###################################
	inline FormalPowerSeries<tfps>&operator*=(const tfps &c){scale(*this,c);return *this;}
	inline FormalPowerSeries<tfps> operator* (const tfps &c)const{FormalPowerSeries<tfps> tmp=*this;return (tmp*=c);}
	// ################################## comparators ###################################
	inline bool operator ==(const FormalPowerSeries<tfps>&B)const{
		bool res=1;rep(i,0,max<int>(size(),B.size()))res=res&&(*this)[i]==B[i];
		return res;
	}inline bool operator !=(const FormalPowerSeries<tfps>&B)const{return !((*this)==B);}
	// ####### shifting operators (x^{shf}*F-> F>>shf or F<<shf depending on shfs sign) #########
	inline FormalPowerSeries<tfps> operator<< (int shf)const{return shift(*this,shf);}
	inline FormalPowerSeries<tfps>&operator<<=(int shf){*this=(*this)<<shf;return *this;}
	inline FormalPowerSeries<tfps> operator>> (int shf)const{return shift(*this,-shf);}
	inline FormalPowerSeries<tfps>&operator>>=(int shf){*this=(*this)>>shf;return *this;}
	// ######## square - derivative - integration - inverse - log - exp - sqrt - euc_div ########
	FormalPowerSeries<tfps>& square(){
		const int Limit=20; if(size()<Limit)
			(*this)=mult_naive(*this,*this);
		else { const int nm=size()*2,lgi=ilog2(nm-1)+1;
			if constexpr(tfps::arbitrary_ntt())
				internal::mult_arbitrary<tfps,bool(1)>(F);
			else{ F.resize(1<<lgi,0);
				internal::fft(F,0);
				for(tfps&ac:F)ac*=ac;
				internal::fft(F,1);
			} trunc(nm);
		} return *this;
	} FormalPowerSeries<tfps> deriv()const{
		FormalPowerSeries<tfps> G(max<int>(1,F.size()-1),0);// G=D(F)
		rep(i,1,F.size())G[i-1]=F[i]*tfps(i);
		return G;
	} FormalPowerSeries<tfps> integ()const{
		FormalPowerSeries<tfps> G(F.size()+1);//D(G)=F
		rep(i,0,F.size())G[i+1]=F[i]/tfps(i+1);
		return G;
	} FormalPowerSeries<tfps> inv(int n) const {
		assert(F[0]);// G*F=1
		FormalPowerSeries<tfps> G={F[0].inv()};//mpow(F[0],mod-2,mod)
		for(int e=2;e<2*n;e<<=1){
			FormalPowerSeries<tfps> ac=G;
			ac*= FormalPowerSeries({F.begin(),F.begin()+min<int>(F.size(),e)}); // gives a ~/2 speedup
			for(tfps &act:ac)act=-act;
			ac+={2}; G*=ac;
			G.trunc(e);
		} return G.trunc(n);
	} FormalPowerSeries<tfps> log(int n) const {
		assert(F[0]==1); //first n terms of G=ln(F)=integ(D(F)/F)
		return (inv(n)*deriv()).trunc(n-1).integ();
	} FormalPowerSeries<tfps> exp(int n) const {
		assert(!F[0]);// first n terms of G=exp(F) 
		FormalPowerSeries<tfps> G={1},ac;
		for(int e=2;e<n*2;e<<=1){
			ac=G.log(e);ac.F.resize(max<int>(ac.size(),min<int>(e,F.size())),0);
			
			rep(i,0,ac.size())ac[i]=(i<int(F.size())?F[i]:tfps(0))- ac[i];
			ac+={1}; (G*=ac).trunc(e);
		} return G.trunc(n);
	} FormalPowerSeries<tfps> pow(lint e,int n)const{
		// first n terms of G=F^e O(nlgn)
		assert(e>=0); if(!e)return {1};
		int xi=0;while(xi<int(F.size())&&F[xi]==0)++xi;
		if( xi>=int(F.size()) || lint(xi) > lint(n)/e )return {0};
		
		// alp could go wrong if its not a integer modulo smthng
		tfps alp=F[xi].pow(e%(tfps::mod()-1)),ainv=F[xi].inv();
		
		int rx=xi?xi*int(e):0; // F^i=exp(i*log(F)) ; shifts and alp,ainv  for making F[0]=1
		return ( ( ( ((*this)>>xi)*ainv).trunc(n-xi).log(n-rx)*tfps(e) 
				).exp(n-rx) >> rx )*alp;
	} FormalPowerSeries<tfps> sqrt(int n)const {
		// this impl only works on numbers modulo a prime
		// first n terms of G^2=F O(nlgn)
		const tfps i2 = tfps(2).inv();//mpow(2,mod-2,mod);
		
		int xi=0;while(xi<int(F.size())&&F[xi]==0)++xi;
		if( xi>=int(F.size()) )return {0};
		
		if(xi&1)return {};
		int vl=internal::Tonelli_Shanks<tfps>(F[xi]);if(vl==-1)return {};
		FormalPowerSeries<tfps> G={vl},ac;//sqrt{F[0]}
		//I assume Tonelli_Shanks returns -1 when it's impossible
		
		FormalPowerSeries<tfps> H=(*this)>>xi;H.trunc(n-xi);
		
		for(int e=2;e<(n-xi)*2;e<<=1){
			ac=G.inv(e)*H;ac.trunc(e);
			G+=ac; G*=i2; G.trunc(e);
		} return (G<<(xi/2)).trunc(n);
	} FormalPowerSeries<tfps> euc_div(const FormalPowerSeries<tfps> &B)const{
		// Given F,B it computes D of -> F=B*D+R
		// Where deg(R)<deg(B) ; note it leaves R easy to compute R=F%B
		FormalPowerSeries<tfps> tf=*this,tb=B;
		tf.trunc(tf.size());tb.trunc(tb.size());
		const int n=tf.size(),m=tb.size(),d=n-m;
		if(d<0)return {1};
		
		auto rev=[&](FormalPowerSeries<tfps>&fn)->FormalPowerSeries<tfps>{
			rep(i,0,fn.size()/2)swap(fn[i],fn[fn.size()-i-1]);
			return fn;
		}; rev(tf); rev(tb);
		tf.trunc(d+1);tb.trunc(d+1);
		return rev((tf*tb.inv(d+1)).trunc(d+1,0)).trunc(d+1);
	}
};typedef FormalPowerSeries<mint> fps;

vector<int> gen_rand(int n){
	vector<int> rs(n);rep(i,0,n){
		rs[i]=rng_64()%mod;
	}
	return rs;
}

int main(){
	ios_base::sync_with_stdio(false);
    cin.tie(NULL);
	cout << setprecision(12) << fixed;
	
	auto eq=[&](const fps &A, const control::fps &B)->bool{
		const int n=max<int>(A.size(),B.size());
		bool res=1;rep(i,0,n){
			int va=A[i],vb=(int(B.size())>i)?int(B[i]):0;
			res&=va==vb;
		}
		return res;
	};
	// t->repetitions;n->coeficients given;k->coefficients generated
	int t,n,k;
	// cin>>t>>n>>k;
	t=100,n=1e5,k=3*1e5;
	// t=1,n=10,k=10;
	
	assert(t>0&&k>0&&n>0);
	rep(sd,0,t){
		vector<int> gen=gen_rand(n);
		// rep(i,0,n)cin>>gen[i];
		fps A;{
			vector<mint> a1(n);rep(i,0,n)a1[i]=gen[i];
			A=fps(a1);}
		control::fps B;{ vector<control::mint>b1(n);
			rep(i,0,n)b1[i]=gen[i];
			B=control::fps(b1);}
		
		// idebug(gen);
		
		if(gen[0]==0)A[0]=1,B[0]=1;
		fps Ai=A.inv(k); control::fps Bi=control::p_inv(B,k);
		
		A[0]=1;B[0]=1;
		fps Al=A.log(k); control::fps Bl=control::p_log(B,k);
		
		A[0]=0;B[0]=0;
		fps Ae=A.exp(k); control::fps Be=control::p_exp(B,k);
		
		A[0]=gen[0];B[0]=gen[0];
		
		lint ex=rng_64();if(ex<0)ex=-ex;
		fps Ap=A.pow(ex,k); control::fps Bp=control::p_pow(B,ex,k);
		
		fps Ap2=A;Ap2.square(); control::fps Bp2=B;p_square(Bp2);
		
		int shf=ex%n;
		fps As1=A<<shf; control::fps Bs1=control::p_shift(B,+shf);
		fps As2=A>>shf; control::fps Bs2=control::p_shift(B,-shf);
		
		fps Asq=A.sqrt(k); control::fps Bsq=control::p_sqrt(B,k);
		Asq.trunc(n);control::p_trunc(Bsq,n);//esta raro el otro caso
		
		fps Ade=A.deriv(); control::fps Bde=control::p_deriv(B);
		fps Ain=A.integ(); control::fps Bin=control::p_inte (B);
		
		bool i=eq(Ai,Bi),l=eq(Al,Bl),e=eq(Ae,Be),p=eq(Ap,Bp),sq=eq(Asq,Bsq),
		s1=eq(As1,Bs1),s2=eq(As2,Bs2),de=eq(Ade,Bde),in=eq(Ain,Bin),
		p2=eq(Ap2,Bp2)
		;
		// idebug(Asq.F);idebug(Bsq);
		// control::p_trunc(Bsq,k);
		// debug(sz(Asq),sz(Bsq));
		// idebug(Asq.F);idebug(Bsq);
		// Asq.square().trunc(n);control::p_square(Bsq);control::p_trunc(Bsq,n);
		// idebug(Asq.F);idebug(A.F);idebug(Bsq);
		// debug(Asq==A,Bsq==B);
		if(!i||!l||!e||!p||!sq||!s1||!s2||!de||!in||!p2){
			
			debug(i,l,e,p,sq,s1,s2,de,in,p2);
			debug(n,k,ex,shf);
			debug(rng_seed);
			assert(0);
		}
		debug("goooooooooooooood",sd);
	}
	return 0;
}