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
const int mod = 998244353;
/*
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
template < uint m, int m_pow=32 >
struct montgomery_int{
	static_assert(m_pow>0&&m_pow<=32); static_assert(m>0);
	static constexpr uint mr=constexpr_calc_mr(m,m_pow);
	constexpr static uint mod(){return m;}
	
	constexpr static uint reduce(ulint x)noexcept{
		
		uint q= uint(x)*mr;
		uint y=((ulint)q*m)>>m_pow;
		
		uint res=(x>>m_pow)+m-y;
		if (res >= m) res -= m;
		return res;
	}
	constexpr static uint transform (uint x)noexcept{return (ulint(x)<<m_pow)%m;}
	constexpr static uint itransform(uint x)noexcept{return reduce(x);}
	
	uint vl;
	
	constexpr montgomery_int()noexcept:vl(0){};
	constexpr montgomery_int( int v)noexcept:vl(transform(v>=0?(v<int(m)?v:v%m):(m+v<m?m+v:(v%m)+m))){};
	constexpr montgomery_int(lint v)noexcept:vl(transform(v>=0?(v<lint(m)?v:v%m):(m-v<m?m-v:(v%m)+m))){};
	constexpr montgomery_int( uint v)noexcept:vl(transform(v<m?v:v%m)){};
	constexpr montgomery_int(ulint v)noexcept:vl(transform(v<ulint(m)?v:v%m)){};
	
	montgomery_int &operator +=(const montgomery_int &ot){ vl+=ot.vl; if(vl>=m)vl-=m; return *this; }
	montgomery_int  operator + (const montgomery_int &ot)const{ return montgomery_int(*this)+=ot; }
	montgomery_int &operator -=(const montgomery_int &ot){ vl=(vl>=ot.vl)?vl-ot.vl:vl+m-ot.vl; return *this; }
	montgomery_int  operator - (const montgomery_int &ot)const{ return montgomery_int(*this)-=ot; }
	montgomery_int &operator *=(const montgomery_int &ot){ vl= reduce((ulint)vl*ot.vl); return *this; }
	montgomery_int  operator * (const montgomery_int &ot)const{ return montgomery_int(*this)*=ot; }
	montgomery_int &operator /=(const montgomery_int &ot){ (*this)*=ot.inv(); return *this; }
	montgomery_int  operator / (const montgomery_int &ot)const{ return montgomery_int(*this)/=ot; }
	
	montgomery_int operator -()const {return montgomery_int(0)-(*this);}
	montgomery_int pow(lint e)const{
		montgomery_int rs=1,ac=*this;while(e){
			if(e&1ll)rs*=ac;
			e>>=1;ac*=ac;
		}
		return rs;
	}
	montgomery_int inv()const{return this->pow(m-2);}//Fermats little theorem
	
	bool operator ==(const montgomery_int &ot)const{return vl==ot.vl;} 
	bool operator !=(const montgomery_int &ot)const{return vl!=ot.vl;}
	bool operator ==(const  uint &ot)const{return reduce(vl)==ot;}
	bool operator ==(const  int &ot )const{return reduce(vl)==uint(ot);}
	
	operator bool() const { return itransform(vl); }
	operator  int() const { return itransform(vl); }
	
	friend ostream &operator<<(ostream &os,const montgomery_int &ac){return os << itransform(ac.vl);}
	friend istream &operator>>(istream&is,montgomery_int &ac){int v;is>>v;ac=montgomery_int(v);return is;}	
}; typedef montgomery_int<mod> mint;


typedef mint tfps; typedef vector<tfps> fps; // fps renaming types

/* START OF NTT */

constexpr int countr_zero_constexpr(unsigned int n) { int x = 0; while (!(n & (1 << x))) x++; return x; }
struct fft_info{
	static const int c_root =3;
	static constexpr int rank2=countr_zero_constexpr(mint::mod()-1);
	array<mint,rank2+1> root,iroot;
	fft_info(){
		root[rank2]=mint(c_root).pow((mint::mod()-1)>>rank2);
		rep(i,rank2-1,-1) root[i]=root[i+1]*root[i+1];
		
		iroot[rank2]=root[rank2].inv();
		rep(i,rank2-1,-1) iroot[i]=iroot[i+1]*iroot[i+1];
	}
};

void fft(fps&A,bool invert){
	const int n=sz(A),lgi=ilog2(n);
	for(int k=0;k<n;++k){
		int j=0,m=k;
		for( int q=0;q<lgi;q+=2	){
			
			tie(j,m)=make_pair(4*j+m%4,m/4);
		}
		if(j>k)debug(j,k);
		if(j>k)swap(A[j],A[k]);
	}
	/*for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1; // counting in reverse
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) swap(A[i], A[j]);
    }*/
	static const fft_info info;
	const array<mint,info.rank2+1> &r_dat=invert?info.iroot:info.root;
	const mint imag=r_dat[2]; 
	debug(imag*imag,imag,mint(-1)+mint(1));
	assert(imag*imag==mint(-1));// i = w_4
	for(int q=2;q<=lgi;q+=2){
		int L=1<<q,r=n/L,Ls=L/4;
		mint w=1,wlen=r_dat[q];
		for(int j=0;j<Ls;++j){
			for(int k=0;k<r;++k){
				mint alpha =      A[k*L     +j];
				mint beta  =w*    A[k*L+  Ls+j];
				mint y     =w*w  *A[k*L+2*Ls+j];
				mint delta =w*w*w*A[k*L+3*Ls+j];
				
				mint r0= alpha + y;
				mint r1= alpha - y;
				mint r2= beta+delta;
				mint r3= beta-delta;
				
				A[k*L     +j] = r0+r2;
				A[k*L+  Ls+j] = r1-imag*r3;//ir3
				A[k*L+2*Ls+j] = r0 - r2;
				A[k*L+3*Ls+j] = r1+imag*r3;// r1+i*r3
				
			}
			w*=wlen;
		}
	}
}

void iter_fft( fps &A,bool invert ){
	// Coley-Tuckey
	const int n=sz(A),lgi=ilog2(n);
	for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1; // counting in reverse
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) swap(A[i], A[j]);
    }
	static const fft_info info;
	// My experimentation shows that (2) scales better
	
	for(int m=1,acb=1;m<=lgi;++m,acb<<=1){
		mint w=1,wlen=(invert)?info.iroot[m]:info.root[m];
		for(int j=0;j<acb;++j) {
			for(int s=0;s<(1<<lgi);s+=1<<m){
				mint u=A[s+j],v=A[s+j+acb]*w;
				A[s+j    ]=u+v;
				A[s+j+acb]=u-v;
			}
			w*=wlen;
		}
	}
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
	int lgi=ilog2(nm-1)+1;if(lgi&1)++lgi;
	A.resize(1<<lgi,0);B.resize(1<<lgi,0);
	
	fps AA=A;
	idebug(A);idebug(AA);
	fft(A,0);//return;
	iter_fft(AA,0);
	idebug(A);idebug(AA);
	// return;
	
	iter_fft(B,0);
	A=AA;
	rep(i,0,sz(A))A[i]*=B[i];
	AA=A;
	
	idebug(A);idebug(AA);
	fft(A,1);
	iter_fft(AA,1);
	idebug(AA);idebug(A);
	
	
	p_trunc(A,nm);
}
void p_mult(fps &A, const fps &B){
	/*auto cnt_zeros=[](const fps &F){int cnt=0;for(tfps ac:F){cnt+=ac==0;if(cnt>=Limit)return false;}return true;};
	else if(cnt_zeros(A)||cnt_zeros(B)) A=p_mult_sparse(A,B);*/
	static const int Limit=0;
	// A=p_mult_naive(A,B);
	if(min(sz(A),sz(B))<=Limit)A=p_mult_naive(A,B);
	else p_mult_fft(A,B);
}
auto take_time=[&](){return std::chrono::high_resolution_clock::now();};
auto get_durat=[&](auto start){ return std::chrono::duration_cast<std::chrono::milliseconds>(take_time() - start).count(); };
int main(){
	ios_base::sync_with_stdio(false);
    cin.tie(NULL);
	cout << setprecision(12) << fixed;
	
	const bool ios=0;
	int n;
	if(ios)cin>>n;
	else scanf("%d",&n);
	
	// lint pm=0;cin>>pm;
	
	fps A(n);
	if(ios)rep(i,0,n)cin>>A[i];
	else rep(i,0,n){int tmp;scanf("%d",&tmp);A[i]=tmp;}
	
	fps B(n);
	if(ios)rep(i,0,n)cin>>B[i];
	else rep(i,0,n){int tmp;scanf("%d",&tmp);B[i]=tmp;}
	
	p_trunc(A,n); p_trunc(B,n);
	
	auto start=take_time();
	p_mult(B,A);
	auto milisec=get_durat(start);
	cerr << "best2radix time " << milisec << '\n';
	if(ios)rep(i,0,2*n)cout << (i<sz(B)?int(B[i]):0) << ' ';
	else rep(i,0,2*n)printf("%d ",(i<sz(B)?int(B[i]):0));
	cout << '\n';
	
	return 0;
}