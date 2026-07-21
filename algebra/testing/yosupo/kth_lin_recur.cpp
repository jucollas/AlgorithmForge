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

mt19937_64 rng_64( chrono::steady_clock::now().time_since_epoch().count() );
constexpr int ilog2( int num ) { return 8*sizeof(int) - __builtin_clz( num ) - 1; }
template<typename tpow,typename texp=lint> constexpr tpow mpow(tpow x,unsigned long long e,tpow m){tpow res=1;while(e){if(e&1)res=(texp(res)*x)%m;e>>=1;x=(texp(x)*x)%m;}return res;}
const int mod=998244353;/*
Author: Oscar Vargas Pabon

It is though to work modulo primes, so inv (.inv,/,/=) and pow may not
    work properly otherwise.

I assume from my template :
inv :: int mpow(int x,int e,int m){int res=1;while(e){if(e&1)res=(res*1ll*x)%m;e>>=1;x=(x*1ll*x)%m;}return res;}    

Tested in testing/test_alghelp.cpp and in fft/ntt stuff
*/
template<__uint64_t raw_m,typename tint=__uint32_t,typename tmul=__uint64_t,bool arbi_ntt=0>
struct modulo_int{ constexpr static tint m=raw_m; static_assert(m>0);
    constexpr static tint mod(){return m;}
    constexpr static bool arbitrary_ntt(){return arbi_ntt;}
    
    tint vl;
    inline constexpr modulo_int()noexcept:vl(0){};
    inline constexpr modulo_int(      int v)noexcept:vl(v>=0?(v<m?v:v%m):(v+m>=0?v+m:(m+(v%m))%m)){};
    inline constexpr modulo_int(long long v)noexcept:vl(v>=0?(v<m?v:v%m):(v+m>=0?v+m:(m+(v%m))%m)){};

    inline constexpr modulo_int(unsigned       int v)noexcept:vl(v<m?v:v%m){};
    inline constexpr modulo_int(unsigned long long v)noexcept:vl(v<m?v:v%m){};
    
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
    
    inline constexpr operator               bool() const{return vl;}
    inline constexpr operator                int() const{return vl;}
    inline constexpr operator unsigned       int() const{return vl;}
    inline constexpr operator          long long() const{return vl;}
    inline constexpr operator unsigned long long() const{return vl;}
    
    friend ostream &operator<<(ostream &os,const modulo_int &ac){return os << ac.vl;}
    friend istream &operator>>(istream&is,modulo_int &ac){int v;is>>v;ac=modulo_int(v);return is;}  
}; typedef modulo_int<mod> mint;
/* Author: Oscar Vargas Pabon
NTT made by yours truly under the help of several books and
    the almighty idea behind the recent progress on:  https://codeforces.com/blog/entry/151162
previously, NTT was taken from the atcoder library: https://github.com/atcoder/ac-library/
Operations on FPS based on code by MarcosK, other people and multiple blogs all around the place
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
 
Dependencies I assume
everything                     :: my modulo_int template
mult_fft,mult_arbitrary,square :: int ilog2( int num ) { return 8*sizeof(int) - __builtin_clz( num ) - 1; }
Note that in this specific version, the ntt somehow gets slower with montgomery_space 
    (its safe but slower smhw).
 
Notes: Some undefined behaviour may ocur whenever F[ind1]=F[ind2] where F.size()<=ind2 as
        F[ind1] can be evaluated first and the resize in F[ind2] may 'overwrite' the
        reference previously established
        
Notes: (27ull<<59)+1 is a 64-bit prime (IntroMathComputational-Shoup pg 484)
        . (549755813881ull<<24)+1 is too (https://codeforces.com/blog/entry/75326)
*/
namespace internal {
template<typename tint> constexpr int countr_zero_constexpr(tint n){if(n==tint(0))return -1;int x=0;while(!(n&(tint(1)<<x)))++x;return x;}
template<typename tint> constexpr tint primitive_root(tint m){ if ((27ull<<59)+1 == m ) return 5;
    if( (549755813881ull<<24)+1 == m || 998244353 == m )return 3;
    tint dec[32]={},dind=0,md=m-1; for(tint i=2;i*i<=md;++i)if(md%i==0){
        dec[dind++]=i; while(md%i==0)md/=i;
    } if(md>1)dec[dind++]=md; //tests for g such that $\forall_{p|(md-1)} g^{(md-1)/p}=1(mod md)$
    bool fnd=0; tint pr=1;while(!fnd){ ++pr; fnd=1; 
        for(tint i=0;i<dind&&fnd;++i){ // by properties, this will always end
            if constexpr ( std::is_same_v<tint,unsigned long long> )
                fnd=mpow<tint,__uint128_t>(pr,(m-1)/dec[i],m)!=1;
            else if constexpr ( std::is_same_v<tint,long long> )
                fnd=mpow<tint,__uint128_t>(pr,(m-1)/dec[i],m)!=1;
            else fnd=mpow<tint>(pr,(m-1)/dec[i],m)!=1; }
    } return pr;// std::cerr <<pr << " _ primitive root" << endl;
} template<class tfps> constexpr auto prec_rank_fft_root() {
    constexpr int rank=countr_zero_constexpr(tfps::mod()-1);
    std::array<tfps,rank+1> root; // precompute the roots
    root[rank]=tfps(primitive_root(tfps::mod())).pow((tfps::mod()-1)>>rank);
    for(int i=rank;i>0;--i)root[i-1]=root[i]*root[i];
    return root;
}template<class tfps> inline constexpr auto rank_fft_root = prec_rank_fft_root<tfps>();
template<typename tfps,long long N> inline constexpr tfps fft_root(){
    if constexpr(N<=0ll||countr_zero_constexpr(N)>=int(rank_fft_root<tfps>.size()))return tfps(0);
    return rank_fft_root<tfps>[countr_zero_constexpr(N)];
}template<typename tfps,long long N> inline constexpr tfps fft_iroot(){return fft_root<tfps,N>().inv();}
template<typename tfps,bool inv,int n> void butterfly_rec(tfps*a){
    constexpr int m=n/2; if constexpr(n<=1)return;
    constexpr tfps wlen=inv?fft_iroot<tfps,n>():fft_root<tfps,n>();
    tfps w=1; for(int i=0;i<m;++i){
        const tfps u=a[i],v=a[i+m];
        a[i]=u+v; a[i+m]=(u-v)*w; w*=wlen;
    }butterfly_rec<tfps,inv,m>(a);butterfly_rec<tfps,inv,m>(a+m);
}template<typename tfps,bool inv=bool(0),int n=min<__uint64_t>(1<<30,(tfps::mod()-1)&(1-tfps::mod()))>
void butterfly(std::vector<tfps>&a){
    if constexpr(n<=1)return;
    if (n==int(a.size()))butterfly_rec<tfps,inv,n>(a.data());
    else butterfly<tfps,inv,n/2>(a);
} template<typename tfps,bool inv,int n> void ibutterfly_rec(tfps*a){
    constexpr int m=n/2; if constexpr(n<=1)return;
    constexpr tfps wlen=inv?fft_root<tfps,n>():fft_iroot<tfps,n>();
    ibutterfly_rec<tfps,inv,m>(a); ibutterfly_rec<tfps,inv,m>(a+m);
    tfps w=1;for(int i=0;i<m;++i){
        const tfps u=a[i],v=w*a[i+m];
        a[i]=u+v; a[i+m]=u-v; w*=wlen; }
}template<typename tfps,bool inv=bool(0),int n=min<__uint64_t>(1<<30,(tfps::mod()-1)&(1-tfps::mod()))>
void ibutterfly(std::vector<tfps>&a){
    if constexpr(n<=1)return;
    if (n==int(a.size()))ibutterfly_rec<tfps,inv,n>(a.data());
    else ibutterfly<tfps,inv,n/2>(a);
} template<typename tfps,int N=min<__uint64_t>(1<<30,(tfps::mod()-1)&(1-tfps::mod()))>
    void fft_doubling(tfps*a,int n){ // given a'[0,n/2) I compute a'[0,n)
    if constexpr(N<=1)return; // assuming a[n/2,n)=0,..,0
    constexpr int N2=N/2; if(N/2<n){
        for(int i=0;i<N2;++i)a[i+N2]=a[i];
        ibutterfly_rec<tfps,bool(0),N2>(a+N2);
        constexpr tfps wlen=fft_root<tfps,N>(),iN2=tfps(N2).inv();
        tfps w=iN2;for(int i=N2;i<N;++i)a[i]*=w,w*=wlen;
        butterfly_rec<tfps,bool(0),N2>(a+N2);
    } else fft_doubling<tfps,N2>(a,n);
} template<typename tfps> inline void fft(std::vector<tfps>&a,bool invert){
    if(invert){ internal::ibutterfly<tfps>(a);
        tfps n_1 = tfps(int(a.size())).inv();
        for (tfps&x:a)x*=n_1;
    }else internal::butterfly<tfps>(a);
} template<typename tfps> inline void transposed_fft(std::vector<tfps> &a,bool invert){
    if(invert){ internal::butterfly<tfps,bool(1)>(a);
        tfps n_1 = tfps(int(a.size())).inv();
        for (tfps&x:a)x*=n_1;
    } else internal::ibutterfly<tfps,bool(1)>(a);
} template<typename tfps,const bool square=0>
void mult_arbitrary(std::vector<tfps>&A,const std::vector<tfps>&B={}){
    // Im guiding miself on nyaan's library for this
    // https://nyaannyaan.github.io/library/ntt/arbitrary-ntt.hpp
    using u128 = __uint128_t;
    const int nm=A.size()+B.size()-1,lgi=ilog2(nm-1)+1,n=1<<lgi,bsz=B.size();
    A.resize(n,tfps(0)); // IntroMathComputational-Shoup pg 484 __ they get overflow
    // constexpr int m0 = (1<<30)*3+1, m1 = (1<<28)*13+1, m2 = (1<<27)*29+1;
    constexpr int m0 = 167772161, m1 = 469762049, m2 = 754974721;
    
    std::vector<modulo_int<m0>> A0(n);
    for(int i=0;i<n;++i)A0[i]=(long long)(A[i]);
    internal::fft<modulo_int<m0>>(A0,0);
    
    if constexpr(square) for(int i=0;i<(1<<lgi);++i)A0[i]*=A0[i];
    else{ std::vector<modulo_int<m0>> B0(n);
        for(int i=0;i<n;++i)B0[i]=(i>=bsz)?0:(long long)(B[i]);
        internal::fft<modulo_int<m0>>(B0,0);
        for(int i=0;i<(1<<lgi);++i)A0[i]*=B0[i];    
    }internal::fft<modulo_int<m0>>(A0,1);
    
    if constexpr ( tfps::mod()<m0 && tfps::mod()*1ll*tfps::mod() < m0 ){
        for(int i=0;i<n;++i)A[i]=int(A0[i]);
        return;
    } std::vector<modulo_int<m1>> A1(1<<lgi);
    for(int i=0;i<n;++i)A1[i]=(long long)(A[i]);
    internal::fft<modulo_int<m1>>(A1,0);
    
    if constexpr(square) for(int i=0;i<(1<<lgi);++i)A1[i]*=A1[i];
    else { std::vector<modulo_int<m1>> B1(1<<lgi);
        for(int i=0;i<n;++i)B1[i]=(i>=bsz)?0:(long long)(B[i]);
        internal::fft<modulo_int<m1>>(B1,0);
        for(int i=0;i<(1<<lgi);++i)A1[i]*=B1[i];
    }internal::fft<modulo_int<m1>>(A1,1);
        
    constexpr bool only_2=0,cond_2=tfps::mod()<std::max<int>(m0,m1) && tfps::mod()*1ll*tfps::mod() < m0*1ll*m1;
    if constexpr ( only_2 || cond_2 ){
        constexpr int i01=modulo_int<m0>(m1).inv(),i10=modulo_int<m1>(m0).inv();
        constexpr long long m01=m0*1ll*m1;
        for(int i=0;i<n;++i) A[i]= (long long)( (
                    int(A0[i])*u128( i01*1ll*m1 ) +
                    int(A1[i])*u128( i10*1ll*m0 ) ) %m01 );
        return;
    } std::vector<modulo_int<m2>> A2(1<<lgi);
    for(int i=0;i<n;++i)A2[i]=(long long)(A[i]);
    internal::fft<modulo_int<m2>>(A2,0);
    if constexpr(square) for(int i=0;i<(1<<lgi);++i)A2[i]*=A2[i];
    else { std::vector<modulo_int<m2>> B2(1<<lgi);
        for(int i=0;i<n;++i)B2[i]=(i>=bsz)?0:(long long)(B[i]);
        internal::fft<modulo_int<m2>>(B2,0);
        for(int i=0;i<(1<<lgi);++i)A2[i]*=B2[i];
    } internal::fft<modulo_int<m2>>(A2,1);
    
    constexpr int i12=modulo_int<m0>(m1*1ll*m2).inv(),
                  i02=modulo_int<m1>(m0*1ll*m2).inv(),
                  i01=modulo_int<m2>(m0*1ll*m1).inv();
    constexpr u128 m012=u128(m0)*m1*m2;
    for(int i=0;i<n;++i) A[i]=(long long)( ( (
              int(A0[i])*u128( (m1*1ll*m2*u128(i12))%m012 ) +
              int(A1[i])*u128( (m0*1ll*m2*u128(i02))%m012 ) +
              int(A2[i])*u128( (m0*1ll*m1*u128(i01))%m012 ) )%m012 )%tfps::mod() );
}  } //end internal namespace
template<typename tfps> struct FormalPowerSeries{ std::vector<tfps> F;
    static FormalPowerSeries<tfps> mult_naive(const FormalPowerSeries<tfps>&A,const FormalPowerSeries<tfps>&B){
        FormalPowerSeries<tfps> C(A.size()+B.size()-1,tfps(0));
        if(A.size()>=B.size())for(int i=0;i<int(A.size());++i)for(int j=0;j<int(B.size());++j)C[i+j]+=A[i]*B[j];
        else                  for(int i=0;i<int(B.size());++i)for(int j=0;j<int(A.size());++j)C[i+j]+=B[i]*A[j];
        return C;
    } static void mult_fft(FormalPowerSeries<tfps>&A,FormalPowerSeries<tfps>B){
        const int nm=A.size()+B.size()-1,lgi=ilog2(nm)+1;// A'=A*B in O(nlgn)
        A.F.resize(1<<lgi,tfps(0));B.F.resize(1<<lgi,tfps(0));
        internal::fft<tfps>(A.F,0);internal::fft<tfps>(B.F,0);
        for(int i=0;i<(1<<lgi);++i) A[i]*=B[i];
        internal::fft<tfps>(A.F,1);
    } static void mult(FormalPowerSeries<tfps>&A, const FormalPowerSeries<tfps>&B){
        const int nm=A.F.size()+B.F.size()-1;
        static const int Limit=20;
        if(std::min<int>(A.size(),B.size())<=Limit)A=mult_naive(A,B);
        else { if constexpr( tfps::arbitrary_ntt() )
                internal::mult_arbitrary<tfps>(A.F,B.F);
            else mult_fft(A,B);
        } A.trunc(nm);
    } static void scale(FormalPowerSeries<tfps>&F,tfps vl){
        for(tfps &ac:F)ac*=vl; // F*vl
    } static void add(FormalPowerSeries<tfps>&A,const FormalPowerSeries<tfps>&B,int sgn=1){
        A.F.resize(std::max<int>(A.size(),B.size()),tfps(0));// A+B*sgn
        if(sgn==-1)    for(int i=0;i<std::min<int>(A.size(),B.size());++i)A[i]-=B[i];
        else if(sgn==1)for(int i=0;i<std::min<int>(A.size(),B.size());++i)A[i]+=B[i];
        else assert(0); // I assume $sgn\in\{-1,1\}$
    } static FormalPowerSeries<tfps> shift(const FormalPowerSeries<tfps>&F,int xi){
        FormalPowerSeries<tfps> G(std::max<int>(1,F.size()+xi),tfps(0)); // G=x^{xi}F
        for(int i=0;i<int(F.size());++i)if(i+xi>=0&&i+xi<int(G.size()))G[i+xi]=F[i];
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
    FormalPowerSeries<tfps>& trunc(int n=-1,bool elim_0=1){ if(n==-1)n=F.size();
        F.resize(std::max<int>(1,std::min<int>(F.size(),n)));
        if(elim_0)while(int(F.size())>1&&F.back()==tfps(0))F.pop_back();
        return *this;
    } int size()const{return F.size();} // the &[] may cause some bugs
    tfps &operator[](int ind){F.resize(std::max<int>(size(),ind+1),tfps(0));return F[ind];}
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
        bool res=1;for(int i=0;i<std::max<int>(size(),B.size());++i)res=res&&(*this)[i]==B[i];
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
        else { const int nm=size()*2-1,lgi=ilog2(nm)+1;
            if constexpr(tfps::arbitrary_ntt())
                internal::mult_arbitrary<tfps,bool(1)>(F);
            else{ F.resize(1<<lgi,0); internal::fft(F,0);
                for(tfps&ac:F)ac*=ac;
                internal::fft(F,1);
            } trunc(nm);
        } return *this;
    } FormalPowerSeries<tfps> deriv()const{
        FormalPowerSeries<tfps> G(std::max<int>(1,size()-1),0);// G=D(F)
        for(int i=1;i<size();++i)G[i-1]=F[i]*tfps(i);
        return G;
    } FormalPowerSeries<tfps> integ()const{
        FormalPowerSeries<tfps> G(size()+1);//D(G)=F
        for(int i=0;i<size();++i)G[i+1]=F[i]/tfps(i+1);
        return G;//precomputing inverses has marginal gain~30miliseconds
    } template<int N>void rec_inv(tfps*FF,tfps*G)const{
        if constexpr(N==1)G[0]=F[0].inv();
        if constexpr(N<=1)return;
        rec_inv<N/2>(FF,G); constexpr int N2=N*2;
        const int k=N<int(F.size())?N:F.size(); 
        for(int i=0;i<k;++i)FF[i]=F[i];
        for(int i=k;i<N2;++i)FF[i]=0;
        internal::butterfly_rec<tfps,bool(0),N2>(FF);
        internal::butterfly_rec<tfps,bool(0),N2>(G);
        for(int i=0;i<N2;++i)G[i]=G[i]+G[i]-G[i]*G[i]*FF[i];
        internal::ibutterfly_rec<tfps,bool(0),N2>(G);
        for(int i=0;i<N ;++i)G[i]*=tfps(N2).inv();
        for(int i=N;i<N2;++i)G[i] =0;
    } template<int N=std::min<long long>(1<<25,(tfps::mod()-1)&(1-tfps::mod()))>
    FormalPowerSeries<tfps> inv(int n) const {
        if constexpr(N<1) return FormalPowerSeries<tfps>();
        constexpr int N2=N/2; if(N2<n){ assert(F[0]);
            std::vector<tfps> FF(N*2),GF(N*2,tfps(0));
            rec_inv<N>(FF.data(),GF.data());
            return FormalPowerSeries<tfps>(GF).trunc(n);
        } else return inv<N2>(n);
    } FormalPowerSeries<tfps> log(int n) const {
        assert(F[0]==1); //first n terms of G=ln(F)=integ(D(F)/F)
        return (inv(n)*deriv()).trunc(n-1).integ();
    } FormalPowerSeries<tfps> exp(int n) const {
        assert(!F[0]);// first n terms of G=exp(F) 
        FormalPowerSeries<tfps> G={1},ac;
        for(int e=2;e<n*2;e<<=1){ ac=G.log(e);
            ac.F.resize(std::max<int>(ac.size(),std::min<int>(e,F.size())),tfps(0));
            for(int i=0;i<ac.size();++i)ac[i]=(i<size()?F[i]:tfps(0))-ac[i];
            ac[0]+=1; (G*=ac).trunc(e);
        } return G.trunc(n);
    } FormalPowerSeries<tfps> pow(long long e,int n)const{
        if(!e)return {1}; // first n terms of G=F^e O(nlgn)
        if(e<0)return inv(n).pow(-e,n);
        int xi=0;while(xi<size()&&F[xi]==0)++xi;
        if(xi>=size()||xi>n/e)return {0};// shifts and alp,ainv  
        const tfps alp=F[xi].pow(e),ainv=F[xi].inv();// for making F[0]=1
        const int rx=xi?xi*int(e):0; // F^i=exp(i*log(F))
        return ( ( ( ((*this)>>xi)*ainv).trunc(n-xi).log(n-rx)*tfps(e) 
                ).exp(n-rx) << rx )*alp;
    } FormalPowerSeries<tfps> euc_div(const FormalPowerSeries<tfps> &B)const{
        FormalPowerSeries<tfps> tf=*this,tb=B; //Computes D of F=B*D+R 
        tf.trunc(tf.size());tb.trunc(tb.size()); // Where deg(R)<deg(B)
        const int n=tf.size(),m=tb.size(),d=n-m;
        if(d<0)return {1};//  ; note it leaves R easy to compute R=F%B
        auto rev=[&](FormalPowerSeries<tfps>&fn)->FormalPowerSeries<tfps>&{
            for(int i=0;i<fn.size()/2;++i)std::swap(fn[i],fn[fn.size()-i-1]);
            return fn;
        }; rev(tf).trunc(d+1); rev(tb).trunc(d+1);
        return rev((tf*tb.inv(d+1)).trunc(d+1,0)).trunc(d+1); }
};typedef FormalPowerSeries<mint> fps;
/*
Author:Oscar Vargas Pabon

Taken from:
    * https://codeforces.com/blog/entry/111862
    * "A Simple and Fast Algorithm for Computing the N-th Term of a Linearly Recurrent Sequence"
        Alin Bostan and Ryuhei Mori
Note that other versions as the one computing a 'window' of terms are still not implemented

I assume fps from my impl

Notes: If I want to generate the first n terms then I can also do it by doing
            (P*Q.inv(n)).trunc(n); for F=P/Q the polinomials            
*/
#define FAST_VER
#ifdef FAST_VER
template<typename tfps> tfps bostanMori(std::vector<tfps>P,std::vector<tfps>Q,__uint64_t k){
    const int d=Q.size(),dlgi=ilog2(d-1)+2,d2=1<<dlgi,d22=d2/2;
    P.resize(d2,tfps(0));Q.resize(d2,tfps(0)); // computes |x^k|P/Q in time O(d*lgd*lgk)
    
    vector<tfps>wrot(d22);{ const tfps wlen=internal::rank_fft_root<tfps>[dlgi].inv();
        wrot[0]=1;for(int i=1;i<d22;++i)wrot[i]=wrot[i-1]*wlen;
        for(int i=0;i<d22;++i){
            int x=0;for(int e=0;e<dlgi-1;++e)x|=(i>>e&1)<<(dlgi-2-e);
            if(x<i)swap(wrot[i],wrot[x]);
        }for(int i=0;i<d22;++i)wrot[i]*=tfps(2).inv();
    } internal::fft(Q,0);internal::fft(P,0); while(k){
        for(int i=0;i<d2;++i)P[i]*=Q[i^1];  
        if(k&1ull) for(int i=0;i<d22;++i)P[i]=(P[i*2]-P[i*2+1])*wrot[i];
        else for(int i=0;i<d22;++i)P[i]=(P[i*2]+P[i*2+1])*tfps(2).inv();
        internal::fft_doubling(P.data(),d2);

        for(int i=0;i<d22;++i)Q[i]=Q[i<<1]*Q[i<<1|1];
        internal::fft_doubling(Q.data(),d2);
        k>>=1;
    } internal::fft(Q,1);internal::fft(P,1);
    return P[0]/Q[0];
}
#else
template<typename tfps> tfps bostanMori(std::vector<tfps>P,std::vector<tfps>Q,__uint64_t k){
    // computes |x^k|P/Q in time O(d*lgd*lgk)
    const int d=Q.size(),d2=1<<(ilog2(d-1)+2);
    std::vector<tfps>nQ(d2);
    P.resize(d2,0);Q.resize(d2,0); while(k){ 
        rep(i,0,d2)nQ[i]=(i&1)?-Q[i]:Q[i];

        internal::fft(Q,0);internal::fft(P,0);internal::fft(nQ,0);
        for(int i=0;i<d2;++i)P[i]*=nQ[i],Q[i]*=nQ[i];//saves 3 ffts
        internal::fft(Q,1);internal::fft(P,1);
        
        for(int i=0;i<d;++i)Q[i]=Q[i*2],P[i]=P[i*2+(k&1ll)];
        k>>=1; for(int i=d;i<d2;++i)Q[i]=P[i]=0;
    } return P[0]/Q[0];
}
#endif
mint seq_bostanMori( fps P,fps Q,lint k){
    // P are assumed to be the first |P| terms of the sequence
    // Q is assumed to be the characteristic polynomial of degree d
    // $F_i=\sum_j F_{i-j}q_j$ then $Q=x^d-\sum_{j=0}^{d-1} q_jx^j$
    //reverse(Q.begin(),Q.end());
    // the reverse transforms to representation
    // $Q=1-\sum_{j=0}^{d-1}x^{j+1}q_j$
    const int d=Q.size();
    P=(P.trunc(d)*Q).trunc(d-1);
    return bostanMori<mint>(P.F,Q.F,k); }

void solve() {
	int d;cin>>d;
    lint k;cin>>k;
    fps A(d),C(d+1);
    rep(i,0,d)cin>>A[i];
    rep(i,1,d+1)cin>>C[i];
    C[0]=1;rep(i,1,d+1)C[i]=-C[i];
    cout << seq_bostanMori(A,C,k) << '\n';
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