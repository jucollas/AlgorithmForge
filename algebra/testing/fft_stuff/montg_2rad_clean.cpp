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

typedef unsigned int uint; typedef unsigned long long ulint;
constexpr uint constexpr_calc_modr(uint n,int m_pow){uint nr=1;for(int i=0;i<ilog2(m_pow);++i)nr*=2-n*+nr;return nr;}
struct montgomery_int{
	static const int m_pow=32;static_assert(m_pow>0&&m_pow<=32);
	static const uint mod= 998244353; static_assert(mod>0);
	static constexpr uint modr=constexpr_calc_modr(mod,m_pow);
	
	constexpr static uint reduce(ulint x)noexcept{
		
		uint q= uint(x)*modr;
		uint m=((ulint)q*mod)>>m_pow;
		
		uint res=(x>>m_pow)+mod-m;
		if (res >= mod) res -= mod;
		return res;
	}
	constexpr static uint transform (uint x)noexcept{return (ulint(x)<<m_pow)%mod;}
	constexpr static uint itransform(uint x)noexcept{return reduce(x);}
	
	uint vl;
	
	constexpr montgomery_int()noexcept:vl(0){};
	constexpr montgomery_int( int v)noexcept:vl(transform(v>=0?v%mod:(v%mod)+mod)){};
	constexpr montgomery_int(lint v)noexcept:vl(transform(v>=0?v%mod:(v%mod)+mod)){};
	constexpr montgomery_int(uint v)noexcept:vl(transform(v%mod)){};
	constexpr montgomery_int(ulint v)noexcept:vl(transform(v%mod)){};
	
	montgomery_int &operator +=(const montgomery_int &ot){ vl+=ot.vl; if(vl>=mod)vl-=mod; return *this; }
	montgomery_int  operator + (const montgomery_int &ot)const{ return montgomery_int(*this)+=ot; }
	montgomery_int &operator -=(const montgomery_int &ot){ vl=(vl>=ot.vl)?vl-ot.vl:vl+mod-ot.vl; return *this; }
	montgomery_int  operator - (const montgomery_int &ot)const{ return montgomery_int(*this)-=ot; }
	montgomery_int &operator *=(const montgomery_int &ot){ vl= reduce((ulint)vl*ot.vl); return *this; }
	montgomery_int  operator * (const montgomery_int &ot)const{ return montgomery_int(*this)*=ot; }
	montgomery_int &operator /=(const montgomery_int &ot){ (*this)*=ot.inverse(); return *this; }
	montgomery_int  operator / (const montgomery_int &ot)const{ return montgomery_int(*this)/=ot; }
	
	montgomery_int operator -()const {return montgomery_int(0)-(*this);}
	montgomery_int pow(lint e)const{
		montgomery_int rs=1,ac=*this;while(e){
			if(e&1ll)rs*=ac;
			e>>=1;ac*=ac;
		}
		return rs;
	}
	montgomery_int inverse()const{return this->pow(mod-2);}//Fermats little theorem
	
	bool operator ==(const montgomery_int &ot)const{return vl==ot.vl;} 
	bool operator ==(const  uint &ot)const{return vl==transform(ot%mod);}
	bool operator ==(const  int &ot )const{return vl==transform(ot%mod);}
	bool operator !=(const montgomery_int &ot)const{return vl!=ot.vl;}
	
	operator bool() const { return itransform(vl); }
	operator  int() const { return itransform(vl); }
	
	friend ostream &operator<<(ostream &os,const montgomery_int &ac){return os << itransform(ac.vl);}
	friend istream &operator>>(istream&is,montgomery_int &ac){int v;cin>>v;ac=montgomery_int(v);return is;}	
}; const int mod=montgomery_int::mod; typedef montgomery_int mint;

typedef mint tfps; typedef vector<tfps> fps; // fps renaming types

/* START OF NTT */

constexpr int countr_zero_constexpr(unsigned int n) { int x = 0; while (!(n & (1 << x))) x++; return x; }
struct fft_info{
	static const int c_root =3;
	static constexpr int rank2=countr_zero_constexpr(mint::mod-1);
	array<mint,rank2+1> root,iroot;
	fft_info(){
		root[rank2]=mint(c_root).pow((mint::mod-1)>>rank2);
		rep(i,rank2-1,-1) root[i]=root[i+1]*root[i+1];
		
		iroot[rank2]=root[rank2].inverse();
		rep(i,rank2-1,-1) iroot[i]=iroot[i+1]*iroot[i+1];
	}
};

void fft(fps&A,bool invert){
	const int n=sz(A),lgi=ilog2(n);
	static const fft_info info;
	if(invert){		
		// invert with coley-tukey
		for(int e=1,acb=1;e<=lgi;++e,acb<<=1){
			mint w=1,wlen=info.iroot[e];
			for(int s=0;s<(1<<lgi);s+=1<<e){
				for(int j=0;j<acb;++j){
					mint u=A[s+j],v=A[s+j+acb]*w;
					A[s+j    ]=u+v;
					A[s+j+acb]=u-v;
					
					w*=wlen;
				}
				w=-w;
			}
		}
		mint n_1 = mint(n).inverse();
        for (mint & x : A)x*=n_1;
	}else{
		// direct with gentlema-sandle
		for(int e=lgi,abs=1<<(e-1);e;--e,abs>>=1){
			mint w=1,wlen=info.root[e];
			for(int j=0;j<(1<<(e-1));++j){
				for(int s=0;s<n;s+=1<<e){
					mint u=A[s+j],v=A[s+j+abs];
					A[s+j]=u+v;
					A[s+j+abs]=(u-v)*w;
				}
				w*=wlen;
			}
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
	int lgi=ilog2(nm-1)+1;A.resize(1<<lgi,0);B.resize(1<<lgi,0);
	fft(A,0);//return;
	fft(B,0);
	rep(i,0,sz(A))A[i]*=B[i];
	fft(A,1);
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