/*
Author: Oscar Vargas Pabon

NTT taken from the atcoder library
	https://github.com/atcoder/ac-library/

Based on code by MarcosK, other people and multiple blogs all around the place
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
Note that in this specific version, the ntt wont work well if using
	montgomery_space. Ill eventually get a montgomery-safe impl

Remember 'primitive_root(mod)' will return a primitive root of the modulo in time
	$O(\sqrt{mod}+R log^2 mod)$ usefull to get the constant

REMEMBER TO COMMENT 'SHITTY_GCC_VERSION' when sending/using arbitrary_modulus multiplication
	as my shitty gcc version doesnt support 128-bit integers or 'if constexpr'

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
Notes: I havent fully finished mult_arbitrary to handle arbitrary modulus.
		I'm missing cases where I use 3 modulus to do the calculations
		
Notes: (1<<59)*27+1 is a 64-bit prime (IntroMathComputational-Shoup pg 484)


*/

/* START OF NTT */
namespace internal { // taken from atcoder

#define SHITTY_GCC_VERSION

int countr_zero(unsigned int n) { return __builtin_ctz(int32_t(n)); }
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
}
template<typename tfps>
struct fft_info {
	static const int g =(tfps::mod()==998244353)?3:primitive_root(int(tfps::mod())); // primitive root 
	
    static constexpr int rank2 = countr_zero_constexpr(int(tfps::mod()) - 1);
    std::array<tfps, rank2 + 1> root;   // root[i]^(2^i) == 1
    std::array<tfps, rank2 + 1> iroot;  // root[i] * iroot[i] == 1

    std::array<tfps, std::max(0, rank2 - 2 + 1)> rate2;
    std::array<tfps, std::max(0, rank2 - 2 + 1)> irate2;

    std::array<tfps, std::max(0, rank2 - 3 + 1)> rate3;
    std::array<tfps, std::max(0, rank2 - 3 + 1)> irate3;

    fft_info() {
        root[rank2] = tfps(g).pow((tfps::mod() - 1) >> rank2);
        iroot[rank2] = root[rank2].inv();
        for (int i = rank2 - 1; i >= 0; i--) {
            root[i] = root[i + 1] * root[i + 1];
            iroot[i] = iroot[i + 1] * iroot[i + 1];
        }

        {
            tfps prod = 1, iprod = 1;
            for (int i = 0; i <= rank2 - 2; i++) {
                rate2[i] = root[i + 2] * prod;
                irate2[i] = iroot[i + 2] * iprod;
                prod *= iroot[i + 2];
                iprod *= root[i + 2];
            }
        }
        {
            tfps prod = 1, iprod = 1;
            for (int i = 0; i <= rank2 - 3; i++) {
                rate3[i] = root[i + 3] * prod;
                irate3[i] = iroot[i + 3] * iprod;
                prod *= iroot[i + 3];
                iprod *= root[i + 3];
            }
        }
    }
};

template<typename tfps>
void butterfly(std::vector<tfps>& a) {
    int n = int(a.size());
    int h = countr_zero((unsigned int)n);

    static const fft_info<tfps> info;

    int len = 0;  // a[i, i+(n>>len), i+2*(n>>len), ..] is transformed
    while (len < h) {
        if (h - len == 1) {
            int p = 1 << (h - len - 1);
            tfps rot = 1;
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
            tfps rot = 1, imag = info.root[2];
            for (int s = 0; s < (1 << len); s++) {
                tfps rot2 = rot * rot;
                tfps rot3 = rot2 * rot;
                int offset = s << (h - len);
                for (int i = 0; i < p; i++) {
                    auto mod2 = 1ULL * tfps::mod() * tfps::mod();
                    auto a0 = 1ULL * a[i + offset].vl;
                    auto a1 = 1ULL * a[i + offset + p].vl * rot.vl;
                    auto a2 = 1ULL * a[i + offset + 2 * p].vl * rot2.vl;
                    auto a3 = 1ULL * a[i + offset + 3 * p].vl * rot3.vl;
                    auto a1na3imag =
                        1ULL * tfps(a1 + mod2 - a3).vl * imag.vl;
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
template<typename tfps>
void butterfly_inv(std::vector<tfps>& a) {
    int n = int(a.size());
    int h = countr_zero((unsigned int)n);

    static const fft_info<tfps> info;

    int len = h;  // a[i, i+(n>>len), i+2*(n>>len), ..] is transformed
    while (len) {
        if (len == 1) {
            int p = 1 << (h - len);
            tfps irot = 1;
            for (int s = 0; s < (1 << (len - 1)); s++) {
                int offset = s << (h - len + 1);
                for (int i = 0; i < p; i++) {
                    auto l = a[i + offset];
                    auto r = a[i + offset + p];
                    a[i + offset] = l + r;
                    a[i + offset + p] =
                        (unsigned long long)((unsigned int)(l.vl - r.vl) + tfps::mod()) *
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
            tfps irot = 1, iimag = info.iroot[2];
            for (int s = 0; s < (1 << (len - 2)); s++) {
                tfps irot2 = irot * irot;
                tfps irot3 = irot2 * irot;
                int offset = s << (h - len + 2);
                for (int i = 0; i < p; i++) {
                    auto a0 = 1ULL * a[i + offset + 0 * p].vl;
                    auto a1 = 1ULL * a[i + offset + 1 * p].vl;
                    auto a2 = 1ULL * a[i + offset + 2 * p].vl;
                    auto a3 = 1ULL * a[i + offset + 3 * p].vl;

                    auto a2na3iimag =
                        1ULL *
                        tfps((tfps::mod() + a2 - a3) * iimag.vl).vl;

                    a[i + offset] = a0 + a1 + a2 + a3;
                    a[i + offset + 1 * p] =
                        (a0 + (tfps::mod() - a1) + a2na3iimag) * irot.vl;
                    a[i + offset + 2 * p] =
                        (a0 + a1 + (tfps::mod() - a2) + (tfps::mod() - a3)) *
                        irot2.vl;
                    a[i + offset + 3 * p] =
                        (a0 + (tfps::mod() - a1) + (tfps::mod() - a2na3iimag)) *
                        irot3.vl;
                }
                if (s + 1 != (1 << (len - 2)))
                    irot *= info.irate3[countr_zero(~(unsigned int)(s))];
            }
            len -= 2;
        }
    }
}
template<typename tfps>
void fft(std::vector<tfps> &A,bool invert){
	if(invert){
		internal::butterfly_inv<tfps>(A);
		tfps n_1 = tfps(sz(A)).inv();
		for (tfps & x : A)x*=n_1;
	} else internal::butterfly<tfps>(A);
}
template<typename tfps>
void transposed_fft(std::vector<tfps> &A,bool invert){
	if (!invert) internal::butterfly_inv<tfps>(A);
	reverse(A.begin() + 1, A.end());
	if(invert){
		internal::butterfly<tfps>(A);
		tfps n_1=tfps(sz(A)).inv();
		for (tfps &x: A) x *= n_1;
	}
}
template<typename tfps>
void mult_arbitrary(std::vector<tfps>&A,const std::vector<tfps> &B ){
	// Im guiding miself on nyaan's library for this
	// https://nyaannyaan.github.io/library/ntt/arbitrary-ntt.hpp
	#ifndef SHITTY_GCC_VERSION
		using u128 = __uint128_t;
	#else 
		using u128=unsigned long long;
	#endif
	const int nm=A.size()+B.size(),lgi=ilog2(nm-1)+1,n=1<<lgi,bsz=B.size();
	A.resize(n,0); // IntroMathComputational-Shoup pg 484 __ they get overflow
	// constexpr int m0 = (1<<30)*3+1, m1 = (1<<28)*13+1, m2 = (1<<27)*29+1;
	constexpr int m0 = 167772161, m1 = 469762049, m2 = 754974721;
	
	std::vector<modulo_int<m0>> A0(n),B0=A0;
	for(int i=0;i<n;++i)A0[i]=lint(A[i]),B0[i]=(i>=bsz)?0:lint(B[i]);
	internal::fft<modulo_int<m0>>(A0,0);internal::fft<modulo_int<m0>>(B0,0);
	rep(i,0,1<<lgi)A0[i]*=B0[i];
	internal::fft<modulo_int<m0>>(A0,1);
	
	#ifndef SHITTY_GCC_VERSION
	if constexpr ( tfps::mod()<m0 && tfps::mod()*1ll*tfps::mod() < m0 ){
	#else
	if ( tfps::mod()<m0 && tfps::mod()*1ll*tfps::mod() < m0 ){
	#endif
		for(int i=0;i<n;++i)A[i]=int(A0[i]);
		return;
	}
	std::vector<modulo_int<m1>> A1(1<<lgi),B1=A1;
	for(int i=0;i<n;++i)A1[i]=lint(A[i]),B1[i]=(i>=bsz)?0:lint(B[i]);
	internal::fft<modulo_int<m1>>(A1,0);internal::fft<modulo_int<m1>>(B1,0);
	rep(i,0,1<<lgi)A1[i]*=B1[i];
	internal::fft<modulo_int<m1>>(A1,1);
		
		
	constexpr bool only_2=1,cond_2=tfps::mod()<max<int>(m0,m1) && tfps::mod()*1ll*tfps::mod() < m0*1ll*m1;
	#ifndef SHITTY_GCC_VERSION
	if constexpr ( only_2 || cond_2 ){
	#else
	if ( only_2 || cond_2 ){
	#endif
		constexpr int i01=modulo_int<m0>(m1).inv(),i10=modulo_int<m1>(m0).inv();
		constexpr lint m01=m0*1ll*m1;
		for(int i=0;i<n;++i){
			A[i]= lint( (
			        int(A0[i])*u128( i01*1ll*m1 ) +
				    int(A1[i])*u128( i10*1ll*m0 ) ) %m01 );
		} return;
	}
	assert(0); // esta parte aun no esta terminada xdxdxdxdxdxd
	
	std::vector<modulo_int<m2>> A2(1<<lgi),B2=A2;
	for(int i=0;i<n;++i)A2[i]=lint(A[i]),B2[i]=(i>=bsz)?0:lint(B[i]);
	internal::fft<modulo_int<m2>>(A2,0);internal::fft<modulo_int<m2>>(B2,0);
	rep(i,0,1<<lgi)A2[i]*=B2[i];
	internal::fft<modulo_int<m2>>(A2,1);
		
	constexpr int i12=modulo_int<m0>(m1*1ll*m2).inv(),
				  i02=modulo_int<m1>(m0*1ll*m2).inv(),
				  i01=modulo_int<m2>(m0*1ll*m1).inv();
	constexpr u128 m012=u128(m0)*m1*m2;
	for(int i=0;i<n;++i){
		A[i]=lint( (
			  int(A0[i])*u128( (m1*1ll*m2*u128(i12))%m012 ) +
			  int(A1[i])*u128( (m0*1ll*m2*u128(i02))%m012 ) +
			  int(A2[i])*u128( (m0*1ll*m1*u128(i01))%m012 ) )%tfps::mod() );
	}
}


template<typename tfps>
int Tonelli_Shanks(tfps a) {
	// usado por sqrt cuando trabajo con enteros modulo algo
	//plagiado epicamente de https://judge.yosupo.jp/submission/270105
	const int mod=tfps::mod();
	if (a < 2) return a;
	if (mpow<int>(a, (mod - 1) / 2, mod) != 1) return -1;
	if (mod % 4 == 3) return mpow<int>(a, (mod + 1) / 4, mod);

	tfps b = 3;
	if (mod != 998244353) {
		while (mpow<int>(b, (mod - 1) / 2, mod) == 1) {
			b=tfps(int(rng_64()%(mod-3)) + 2);
		}
	}

	int q = mod - 1,Q = 0;
	while ( !(q&1) ) Q++, q /= 2;

	tfps x = mpow<int>(a, (q + 1) / 2, mod);
	b = mpow<int>(b, q, mod);

	int shift = 2;
	while ( x*x != a) {
		tfps error= tfps(mpow<int>(a,mod-2,mod))*x*x;
		if (mpow<int>(error, 1 << (Q - shift), mod) != 1) x *= b;
		b *=b;
		++shift;
	}
	return x;
}
} //end internal namespace

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
		internal::fft<tfps>(A.F,0);internal::fft<tfps>(B.F,0);
		rep(i,0,sz(A))A[i]*=B[i];
		internal::fft<tfps>(A.F,1);
	}
	
	static void mult(FormalPowerSeries<tfps> &A, const FormalPowerSeries<tfps> &B){
		const int nm=A.F.size()+B.F.size();
		static const int Limit=0;
		if(min(sz(A),sz(B))<=Limit)A=mult_naive(A,B);
		else {
			#ifndef SHITTY_GCC_VERSION
			if constexpr( tfps::arbitrary_ntt() )
			#else
			if ( tfps::arbitrary_ntt() )
			#endif
				internal::mult_arbitrary<tfps>(A.F,B.F);
			else mult_fft(A,B);
		} A.trunc(nm);
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
		return G;
	}
	
	// constructor
	constexpr FormalPowerSeries()noexcept:F({0}){};
	constexpr FormalPowerSeries(const std::vector<tfps> &Fp):F(Fp){};
	FormalPowerSeries(std::initializer_list<tfps> Fp):F(Fp){};
	constexpr FormalPowerSeries(int n,tfps vl=tfps(0)){F.resize(n,vl);}
	
	// iterator stuff
	using iterator=typename std::vector<tfps>::iterator;
	using const_iterator=typename std::vector<tfps>::const_iterator;
	iterator begin() { return F.begin(); } iterator end() { return F.end(); }
	const_iterator begin() const { return F.begin(); } const_iterator end() const { return F.end(); }
	
	// utilities
	FormalPowerSeries<tfps>& trunc(int n=-1,bool elim_0=1){
		if(n==-1)n=sz(F);
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
	
	FormalPowerSeries<tfps> operator % (const FormalPowerSeries<tfps>&B)const{return ((*this)-B*euc_div(B)).trunc(sz(B)); }
	FormalPowerSeries<tfps> &operator %=(const FormalPowerSeries<tfps>&B)const{ (*this)=(*this)%B; return *this;}
	
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
	
	// square - derivative - integration - inverse - log - exp - sqrt - euc_div
	
	FormalPowerSeries<tfps>& square(){
		const int Limit=50;
		if(sz(F)<Limit)(*this)=mult_naive(*this,*this);
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
			
			rep(i,0,sz(ac))ac[i]=(i<sz(F)?F[i]:tfps(0))- ac[i];
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
		tfps alp=F[xi].pow(e%(tfps::mod()-1)),ainv=F[xi].inv();
		
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
		int vl=internal::Tonelli_Shanks<tfps>(F[xi]);if(vl==-1)return {};
		FormalPowerSeries<tfps> G={vl},ac;//sqrt{F[0]}
		//I assume Tonelli_Shanks returns -1 when it's impossible
		
		FormalPowerSeries<tfps> H=(*this)>>xi;H.trunc(n-xi);
		
		for(int e=2;e<(n-xi)*2;e<<=1){
			ac=G.inv(e)*H;ac.trunc(e);
			G+=ac; G*=i2; G.trunc(e);
		}
		return (G<<(xi/2)).trunc(n);
	}
	FormalPowerSeries<tfps> euc_div(const FormalPowerSeries<tfps> &B)const{
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