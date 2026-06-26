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
#include <cassert>
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
template<typename tpow> constexpr tpow mpow(tpow x,lint e,tpow m){tpow res=1;while(e){if(e&1ll)res=(res*1ll*x)%m;e>>=1;x=(x*1ll*x)%m;}return res;}
const int template_limit = 1e6;
// int a[template_limit], b[template_limit];

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
everything     :: #define sz(vec) int(vec.size())
mult,square    :: int ilog2( int num ) { return 8*sizeof(int) - __builtin_clz( num ) - 1; }
TonelliShanks  :: mt19937_64 rng_64( chrono::steady_clock::now().time_since_epoch().count() );
sqrt           :: TonelliShanks, tfps(?).inv()
pow            :: tfps(?).inv(), tfps(?).pow()
 
Notes: Some undefined behaviour may ocur whenever F[ind1]=F[ind2] where sz(F)<=ind2 as
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
        tfps n_1 = tfps(sz(a)).inv();
        for (tfps&x:a)x*=n_1;
	}else internal::butterfly(a);
} template<typename tfps> void transposed_fft(std::vector<tfps> &A,bool invert){
	if (!invert) internal::ibutterfly<tfps>(A);
	reverse(A.begin() + 1, A.end()); if(invert){
		internal::butterfly<tfps>(A);
		tfps n_1=tfps(sz(A)).inv();
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
		FormalPowerSeries<tfps> C(sz(A)+sz(B),0);
		if(sz(A)>=sz(B)) rep(i,0,sz(A))rep(j,0,sz(B)) C[i+j]+=A[i]*B[j];
		else             rep(i,0,sz(B))rep(j,0,sz(A)) C[i+j]+=B[i]*A[j];
		return C;
	} static void mult_fft(FormalPowerSeries<tfps> &A, FormalPowerSeries<tfps> B){
		// A'=A*B in O(nlgn)
		const int nm=A.size()+B.size();
		const int lgi=ilog2(nm-1)+1;
		A.F.resize(1<<lgi,0);B.F.resize(1<<lgi,0);
		internal::fft<tfps>(A.F,0);internal::fft<tfps>(B.F,0);
		rep(i,0,sz(A))A[i]*=B[i];
		internal::fft<tfps>(A.F,1);
	} static void mult(FormalPowerSeries<tfps> &A, const FormalPowerSeries<tfps> &B){
		const int nm=A.F.size()+B.F.size();
		static const int Limit=70;
		if(min(sz(A),sz(B))<=Limit)A=mult_naive(A,B);
		else {
			if constexpr( tfps::arbitrary_ntt() )
				internal::mult_arbitrary<tfps>(A.F,B.F);
			else mult_fft(A,B);
		} A.trunc(nm);
	} static void scale(FormalPowerSeries<tfps> &F,tfps vl){
		for(tfps &ac:F)ac*=vl; // F*vl
	} static void add(FormalPowerSeries<tfps> &A, const FormalPowerSeries<tfps> &B, int sgn=1 ){
		A.F.resize(max(sz(A),sz(B))); // A + B*sgn ; I assume sgn\in\{-1,1\}
		if(sgn==1)      rep(i,0,min(sz(A),sz(B)))A[i]+=B[i];
		else if(sgn==-1)rep(i,0,min(sz(A),sz(B)))A[i]-=B[i];
		else assert(0);
	} static FormalPowerSeries<tfps> shift(const FormalPowerSeries<tfps> &F,int xi){
		FormalPowerSeries<tfps> G(max(1,sz(F)+xi),0); // G=x^{xi}F
		rep(i,0,sz(F))if(i+xi>=0&&i+xi<sz(G)) G[i+xi]=F[i];
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
		if(n==-1)n=sz(F);
		F.resize(max(1,min(sz(F),n)));
		if(elim_0)while(sz(F)>1&&F.back()==0)F.pop_back();
		return *this;
	} int size()const{return F.size();}
	tfps &operator[](int ind){F.resize(max(size(),ind+1),0);return F[ind];}
	tfps operator[](int ind)const {return (size()>ind)?F[ind]:tfps(0);}
	// ############################ basic operators ###############################
	inline FormalPowerSeries<tfps>&operator*=(const FormalPowerSeries<tfps>&B){mult(*this,B); return *this;}
	inline FormalPowerSeries<tfps> operator* (const FormalPowerSeries<tfps>&B)const{FormalPowerSeries<tfps> tmp=*this; return (tmp*=B);}
	inline FormalPowerSeries<tfps>&operator+=(const FormalPowerSeries<tfps>&B){add(*this,B); return *this;}
	inline FormalPowerSeries<tfps> operator+ (const FormalPowerSeries<tfps>&B)const{FormalPowerSeries<tfps> tmp=*this; return (tmp+=B);}
	inline FormalPowerSeries<tfps>&operator-=(const FormalPowerSeries<tfps>&B){add(*this,B,-1); return *this;}
	inline FormalPowerSeries<tfps> operator- (const FormalPowerSeries<tfps>&B)const{FormalPowerSeries<tfps> tmp=*this; return (tmp-=B);}
	inline FormalPowerSeries<tfps> operator% (const FormalPowerSeries<tfps>&B)const{return ((*this)-B*euc_div(B)).trunc(sz(B)); }
	inline FormalPowerSeries<tfps>&operator%=(const FormalPowerSeries<tfps>&B)const{ (*this)=(*this)%B; return *this;}
	inline FormalPowerSeries<tfps> operator-()const{return FormalPowerSeries<tfps>({0}) - (*this);}
	//############################### scalar operators ###################################
	inline FormalPowerSeries<tfps>&operator*=(const tfps &c){scale(*this,c);return *this;}
	inline FormalPowerSeries<tfps> operator* (const tfps &c)const{FormalPowerSeries<tfps> tmp=*this;return (tmp*=c);}
	// ################################## comparators ###################################
	inline bool operator ==(const FormalPowerSeries<tfps>&B)const{
		bool res=1;rep(i,0,max(size(),B.size()))res=res&&(*this)[i]==B[i];
		return res;
	}inline bool operator !=(const FormalPowerSeries<tfps>&B)const{return !((*this)==B);}
	// ####### shifting operators (x^{shf}*F-> F>>shf or F<<shf depending on shfs sign) #########
	inline FormalPowerSeries<tfps> operator<< (int shf)const{return shift(*this,shf);}
	inline FormalPowerSeries<tfps>&operator<<=(int shf){*this=(*this)<<shf;return *this;}
	inline FormalPowerSeries<tfps> operator>> (int shf)const{return shift(*this,-shf);}
	inline FormalPowerSeries<tfps>&operator>>=(int shf){*this=(*this)>>shf;return *this;}
	// ######## square - derivative - integration - inverse - log - exp - sqrt - euc_div ########
	FormalPowerSeries<tfps>& square(){
		const int Limit=70; if(size()<Limit)
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
		FormalPowerSeries<tfps> G(max(1,sz(F)-1),0);// G=D(F)
		rep(i,1,sz(F))G[i-1]=F[i]*tfps(i);
		return G;
	} FormalPowerSeries<tfps> integ()const{
		FormalPowerSeries<tfps> G(sz(F)+1);//D(G)=F
		rep(i,0,sz(F))G[i+1]=F[i]/tfps(i+1);
		return G;
	} FormalPowerSeries<tfps> inv(int n) const {
		assert(F[0]);// G*F=1
		FormalPowerSeries<tfps> G={F[0].inv()};//mpow(F[0],mod-2,mod)
		for(int e=2;e<2*n;e<<=1){
			FormalPowerSeries<tfps> ac=G;
			ac*= FormalPowerSeries({F.begin(),F.begin()+min(sz(F),e)}); // gives a ~/2 speedup
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
			ac=G.log(e);ac.F.resize(max(sz(ac),min(e,sz(F))),0);
			
			rep(i,0,sz(ac))ac[i]=(i<sz(F)?F[i]:tfps(0))- ac[i];
			ac+={1}; (G*=ac).trunc(e);
		} return G.trunc(n);
	} FormalPowerSeries<tfps> pow(lint e,int n)const{
		// first n terms of G=F^e O(nlgn)
		assert(e>=0); if(!e)return {1};
		int xi=0;while(xi<sz(F)&&F[xi]==0)++xi;
		if( xi>=sz(F) || lint(xi) > lint(n)/e )return {0};
		
		// alp could go wrong if its not a integer modulo smthng
		tfps alp=F[xi].pow(e%(tfps::mod()-1)),ainv=F[xi].inv();
		
		int rx=xi?xi*int(e):0; // F^i=exp(i*log(F)) ; shifts and alp,ainv  for making F[0]=1
		return ( ( ( ((*this)>>xi)*ainv).trunc(n-xi).log(n-rx)*tfps(e) 
				).exp(n-rx) >> rx )*alp;
	} FormalPowerSeries<tfps> sqrt(int n)const {
		// this impl only works on numbers modulo a prime
		// first n terms of G^2=F O(nlgn)
		const tfps i2 = tfps(2).inv();//mpow(2,mod-2,mod);
		
		int xi=0;while(xi<sz(F)&&F[xi]==0)++xi;
		if( xi>=sz(F) )return {0};
		
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
		tf.trunc(sz(tf));tb.trunc(sz(tb));
		const int n=sz(tf),m=sz(tb),d=n-m;
		if(d<0)return {1};
		
		auto rev=[&](FormalPowerSeries<tfps>&fn)->FormalPowerSeries<tfps>{
			rep(i,0,sz(fn)/2)swap(fn[i],fn[sz(fn)-i-1]);
			return fn;
		}; rev(tf); rev(tb);
		tf.trunc(d+1);tb.trunc(d+1);
		return rev((tf*tb.inv(d+1)).trunc(d+1,0)).trunc(d+1);
	}
};typedef FormalPowerSeries<mint> fps;
auto take_time=[](){return std::chrono::high_resolution_clock::now();};
auto get_durat=[](auto start){ return std::chrono::duration_cast<std::chrono::nanoseconds>(take_time() - start).count(); };
//std::chrono::milliseconds
//std::chrono::microseconds
//std::chrono::nanoseconds
void multiplication() {
	const int m=1e5;
	// const int m=100;
	const int u=3,d=120;
	rep(i,u,d){
		fps A(i);rep(j,0,i)A[j]=int(rng_64()%mod);
		fps B(m);rep(j,0,m)B[j]=int(rng_64()%mod);
		
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
/*void expon(){
	const int n=1e1,m_pow=10;
	fps F(n);rep(i,0,n)F[i]=int(rng_64()%mod);
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
	
}*/

int32_t main(){
	ios_base::sync_with_stdio(false);
    cin.tie(NULL);
	cout << setprecision(12) << fixed;

    multiplication();
	// expon();
	return 0;
}

