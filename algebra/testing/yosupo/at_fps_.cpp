#include <bits/stdc++.h>

typedef long long lint;

using namespace std;

#define debug(args...) { string _s = #args; replace(_s.begin(), _s.end(), ',', ' '); stringstream _ss(_s); istream_iterator<string> _it(_ss); raw_debug(_it, args);}
void raw_debug(istream_iterator<string> it) {cerr<<endl;}
template<typename T, typename... Args>
void raw_debug(istream_iterator<string> it, T a, Args... args) { cerr <<"<"<< *it << "->" << a << "> "; raw_debug(++it, args...); }
#define idebug(v) {cout<<'['<<#v<<']';for(const auto &el:v)cout << ' ' << el; cout << endl;}
#define adebug(ar,n) {cout<<'['<<#ar<<']';for(int i=0;i<n;++i)cout << ' ' << ar[i]; cout << endl;}

#define rep(i,strt,end) for(int i = strt ; i !=int(end) ; (int(strt)<int(end))?++i:--i )
// #define rall(vec) vec.rbegin(), vec.rend()
// #define all(vec) vec.begin(), vec.end()
#define sz(vec) int(vec.size())
// #define pb push_back
#define pob pop_back
// #define pf push_front
// #define pof pop_front

mt19937_64 rng_64( chrono::steady_clock::now().time_since_epoch().count() );
constexpr int ilog2( int num ) { return 8*sizeof(int) - __builtin_clz( num ) - 1; }
constexpr int mpow(int x,int e,int m){int res=1;while(e){if(e&1)res=(res*1ll*x)%m;e>>=1;x=(x*1ll*x)%m;}return res;}

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

int main(){
	ios_base::sync_with_stdio(false);
    cin.tie(NULL);
	cout << setprecision(12) << fixed;
	int n;cin>>n;
	
	// lint pm=0;cin>>pm;
	
	fps A(n);rep(i,0,n)cin>>A[i];
	p_trunc(A,n);
	
	// int i1=mpow(2,pm%mod,mod),i2=mpow(2,pm%(mod-1ll),mod);
	// debug(i1,i2);
	
	// fps B=p_inv(A,n);
	// fps B=p_log(A,n);
	// fps B=p_exp(A,m,n);
	//fps B=p_pow(A,pm,n);
	// fps C=p_binpow(A,pm,n); idebug(C);
	// fps D=p_binpow(A,pm%(mod-1),n); idebug(D);
	// fps E=p_binpow(A,mod+1,n);idebug(E);
	
	fps B = p_sqrt(A,n);if(B.empty()){cout << "-1\n";return 0;}
	// fps C=B;p_square(C);p_trunc(C,n);idebug(C);
	rep(i,0,n)cout << (i<sz(B)?B[i].vl:0) << ' ';
	cout << '\n';
	
	return 0;
}