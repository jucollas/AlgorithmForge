#include <bits/stdc++.h>
#include<cassert>
typedef long long lint;

using namespace std;

#define debug(args...) { string _s = #args; replace(_s.begin(), _s.end(), ',', ' '); stringstream _ss(_s); istream_iterator<string> _it(_ss); raw_debug(_it, args);}
void raw_debug(istream_iterator<string> it) {cerr<<endl;}
template<typename T, typename... Args>
void raw_debug(istream_iterator<string> it, T a, Args... args) { cerr <<"<"<< *it << "->" << a << "> "; raw_debug(++it, args...); }
#define idebug(v) {cerr<<'['<<#v<<']';for(const auto &el:v)cerr << ' ' << el; cerr << endl;}
#define adebug(ar,n) {cerr<<'['<<#ar<<']';for(int i=0;i<n;++i)cerr << ' ' << ar[i]; cerr << endl;}

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
template<typename tpow> constexpr tpow mpow(tpow x,lint e,tpow m){tpow res=1;while(e){if(e&1ll)res=(res*1ll*x)%m;e>>=1;x=(x*1ll*x)%m;}return res;}

const int mod= 998244353; /*
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
}; /*
Author: Oscar Vargas Pabon

It is though to work modulo primes, so inv (.inv,/,/=) and pow may not
    work properly otherwise.

I assume from my template :
inv :: int mpow(int x,int e,int m){int res=1;while(e){if(e&1)res=(res*1ll*x)%m;e>>=1;x=(x*1ll*x)%m;}return res;}    

Tested in testing/test_alghelp.cpp and in fft/ntt stuff

Multiplications cost less (as we dont take costly modulus), but
    transforming to montgomery space now implies taking a modulus and
    detransforming from montgomery space implies the same work as a 
    multiplication
    
Tener cuidado si tengo que comparar con muchos enteros o cosas por el estilo

Adapted from https://en.algorithmica.org/hpc/number-theory/montgomery/
*/

typedef unsigned int uint; typedef unsigned long long ulint;
constexpr uint constexpr_calc_mr(uint n,int m_pow){uint nr=1;for(int i=0;i<ilog2(m_pow);++i)nr*=2-n*+nr;return nr;}
template < ulint raw_m, int m_pow=32, typename tint=uint,bool arbi_ntt=0 >
struct montgomery_int{
    constexpr static tint m=raw_m; static_assert(m>0);
    static_assert(m_pow>0&&m_pow<=32); 
    static constexpr tint mr=constexpr_calc_mr(m,m_pow);
    inline constexpr static tint mod(){return m;}
    inline constexpr static bool arbitrary_ntt(){return arbi_ntt;}
    
    inline constexpr static tint reduce(ulint x)noexcept{
        
        tint q= tint(x)*mr;
        tint y=((ulint)q*m)>>m_pow;
        
        tint res=(x>>m_pow)+m-y;
        if (res >= m) res -= m;
        return res;
    }
    inline constexpr static tint transform (tint x)noexcept{return (ulint(x)<<m_pow)%m;}
    inline constexpr static tint itransform(tint x)noexcept{return reduce(x);}
    
    tint vl;
    
    inline constexpr montgomery_int()noexcept:vl(0){};
    inline constexpr montgomery_int( int v)noexcept:vl(transform(v>=0?(v< int(m)?v:v%m):(m+v>=0?m+v:(v%m)+m))){};
    inline constexpr montgomery_int(lint v)noexcept:vl(transform(v>=0?(v<lint(m)?v:v%m):(m+v>=0?m+v:(v%m)+m))){};
    inline constexpr montgomery_int( uint v)noexcept:vl(transform(v<m?v:v%m)){};
    inline constexpr montgomery_int(ulint v)noexcept:vl(transform(v<ulint(m)?v:v%m)){};
    
    inline constexpr montgomery_int&operator +=(const montgomery_int &ot){ vl+=ot.vl; if(vl>=m)vl-=m; return *this; }
    inline constexpr montgomery_int operator + (const montgomery_int &ot)const{ return montgomery_int(*this)+=ot; }
    inline constexpr montgomery_int&operator -=(const montgomery_int &ot){ vl=(vl>=ot.vl)?vl-ot.vl:vl+m-ot.vl; return *this; }
    inline constexpr montgomery_int operator - (const montgomery_int &ot)const{ return montgomery_int(*this)-=ot; }
    inline constexpr montgomery_int&operator *=(const montgomery_int &ot){ vl= reduce((ulint)vl*ot.vl); return *this; }
    inline constexpr montgomery_int operator * (const montgomery_int &ot)const{ return montgomery_int(*this)*=ot; }
    inline constexpr montgomery_int&operator /=(const montgomery_int &ot){ (*this)*=ot.inv(); return *this; }
    inline constexpr montgomery_int operator / (const montgomery_int &ot)const{ return montgomery_int(*this)/=ot; }
    
    inline constexpr montgomery_int operator -()const {return montgomery_int(0)-(*this);}
    inline constexpr montgomery_int pow(lint e)const{
        montgomery_int rs=1,ac=*this;while(e){
            if(e&1ll)rs*=ac;
            e>>=1;ac*=ac;
        }
        return rs;
    }
    inline constexpr montgomery_int inv()const{return this->pow(m-2);}//Fermats little theorem
    
    inline constexpr bool operator ==(const montgomery_int &ot)const{return vl==ot.vl;} 
    inline constexpr bool operator !=(const montgomery_int &ot)const{return vl!=ot.vl;}
    inline constexpr bool operator ==(const  uint &ot)const{return reduce(vl)==ot;}
    inline constexpr bool operator ==(const  int &ot )const{return reduce(vl)==ot;}
    
    inline constexpr operator bool() const { return itransform(vl); }
    inline constexpr operator  int() const { return itransform(vl); }
    inline constexpr operator long long() const{ return int(*this); }
    
    friend ostream &operator<<(ostream &os,const montgomery_int &ac){return os << itransform(ac.vl);}
    friend istream &operator>>(istream&is,montgomery_int &ac){int v;is>>v;ac=montgomery_int(v);return is;}  
};

#ifdef MONTGOMERY
    typedef montgomery_int<mod> mint;
#else
    typedef modulo_int<mod> mint;
#endif


typedef mint tfps; typedef vector<tfps> fps; // fps renaming types

/* START OF NTT */
namespace internal { // taken from atcoder

constexpr int countr_zero_constexpr(unsigned int n) { int x = 0; while (!(n & (1 << x))) x++; return x; }


// #define FFT_INFO_VER
#ifdef FFT_INFO_VER
template<typename T>
struct fft_info{
    static constexpr T c_root =3;
    static constexpr int rank2=countr_zero_constexpr(mint::mod()-1);
    array<T,rank2+1> root,iroot;
    fft_info(){
        root[rank2]=c_root.pow((mint::mod()-1)>>rank2);
        rep(i,rank2-1,-1) root[i]=root[i+1]*root[i+1];
        
        iroot[rank2]=root[rank2].inv();
        rep(i,rank2-1,-1) iroot[i]=iroot[i+1]*iroot[i+1];
    }
};
#else
template<typename T,int rank> constexpr T fft_root(){
    constexpr T c_root=3;
    constexpr int rank2=countr_zero_constexpr(T::mod()-1);
    T res; if constexpr(rank==rank2)
        res=c_root.pow((T::mod()-1)>>rank2);
    else res=fft_root<T,rank+1>()*fft_root<T,rank+1>();
    return res;
} template<typename T,int rank> constexpr T fft_iroot(){
    constexpr T c_root=3;
    constexpr int rank2=countr_zero_constexpr(T::mod()-1);
    T res; if constexpr(rank==rank2)
        res=fft_root<T,rank>().inv();
    else res=fft_iroot<T,rank+1>()*fft_iroot<T,rank+1>();
    return res;
}
// this is just for debugging
template<typename T,int rank> constexpr std::vector<T> get_root(){
    constexpr int rank2=countr_zero_constexpr(T::mod()-1);
    if constexpr(rank==rank2)return {fft_root<T,rank>()};
    else{
        std::vector<T> res=get_root<T,rank+1>();
        res.push_back(fft_root<T,rank>());
        return res;
    }
}
#endif//FFT_INFO_VER

template<typename T,int n> void butterfly_rec(T*a){
    if constexpr(n==1)return;
    constexpr int e=ilog2(n),m=n/2;
    
    #ifdef FFT_INFO_VER
    static const T wlen=fft_info<T>().root[e];
    #else
    constexpr T wlen=fft_root<T,e>();
    #endif //FFT_INFO_VER

    // for(int i=0;i<n;++i)cerr << -a[i] << " ";cerr << endl;
    T w=1; for(int i=0;i<m;++i){
        const T u=a[i],v=a[i+m];
        a[i  ]=u+v;
        a[i+m]=(u-v)*w;
        w*=wlen;

        // const T u=a[i],v=w*a[i+m];
        // a[i  ]=u+v;
        // a[i+m]=u-v;
        // w*=wlen;
    }
    butterfly_rec<T,m>(a  );
    butterfly_rec<T,m>(a+m);
}template<typename T,int n=(1<<countr_zero_constexpr(T::mod()-1))>
void butterfly(std::vector<T>&a){
    if constexpr(n==1)return;
    if (n==int(a.size()))butterfly_rec<T,n>(a.data());
    else butterfly<T,n/2>(a);
} template<typename T,int n> void ibutterfly_rec(T*a){
    if constexpr(n==1)return;
    constexpr int e=ilog2(n),m=n/2;
    
    #ifdef FFT_INFO_VER
    static const T wlen=fft_info<T>().iroot[e];
    #else
    constexpr T wlen=fft_iroot<T,e>();
    #endif //FFT_INFO_VER

    ibutterfly_rec<T,m>(a  );
    ibutterfly_rec<T,m>(a+m);
    T w=1;for(int i=0;i<m;++i){
        const T u=a[i],v=w*a[i+m];
        a[i  ]=u+v;
        a[i+m]=u-v;
        w*=wlen;
    }
}template<typename T,int n=(1<<countr_zero_constexpr(T::mod()-1))>
void ibutterfly(std::vector<T>&a){
    if constexpr(n==1)return;
    if (n==int(a.size()))ibutterfly_rec<T,n>(a.data());
    else ibutterfly<T,n/2>(a);
}

} //end internal namespace
template<typename T>
void fft(vector<T>&a,bool invert){
    const int n=a.size();
    // for (int i = 1, j = 0; i < n; i++) {
    //     int bit = n >> 1;
    //     for (; j & bit; bit >>= 1)
    //         j ^= bit;
    //     j ^= bit;

    //     if (i < j)
    //         swap(a[i], a[j]);//,swap(nm[i],nm[j]);
    // } 
	if(invert){ internal::ibutterfly(a);
        T n_1 = T(sz(a)).inv();
        for (T&x:a)x*=n_1;
	}else internal::butterfly(a);
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
	static const int Limit=0;
	if(min(sz(A),sz(B))<=Limit)A=p_mult_naive(A,B);
	else p_mult_fft(A,B);
}
auto take_time=[](){return std::chrono::high_resolution_clock::now();};
auto get_durat=[](auto start){ return std::chrono::duration_cast<std::chrono::milliseconds>(take_time() - start).count(); };
int main(){
	ios_base::sync_with_stdio(false);
    cin.tie(NULL);
	cout << setprecision(12) << fixed;
	
	const bool ios=0;
	int n;
	if(ios)cin>>n;
	else scanf("%d",&n);

    // auto vc=internal::get_root<mint,0>();
    // reverse(vc.begin(),vc.end());
    // for(auto ac:vc)cerr << ac << " ";
    //     cerr << endl;
	
	// lint pm=0;cin>>pm;
	
	fps A(n);
	if(ios)rep(i,0,n)cin>>A[i];
	else rep(i,0,n){int tmp;scanf("%d",&tmp);A[i]=tmp;}

    // fft(A,0);for(auto ac:A)cerr << ac << " ";cerr << endl; return 0;
	
	fps B(n);
	if(ios)rep(i,0,n)cin>>B[i];
	else rep(i,0,n){int tmp;scanf("%d",&tmp);B[i]=tmp;}
	
	p_trunc(A,n); p_trunc(B,n);
	
	auto start=take_time();
	p_mult(B,A);
	auto milisec=get_durat(start);
    #ifdef MONTGOMERY
	   cerr << "weird-recursive-montgomery time " << milisec << '\n';
    #else
       cerr << "weird-recursive time " << milisec << '\n';
    #endif
	if(ios)rep(i,0,2*n)cout << (i<sz(B)?int(B[i]):0) << ' ';
	else rep(i,0,2*n)printf("%d ",(i<sz(B)?int(B[i]):0));
	cout << '\n';
	
	return 0;
}