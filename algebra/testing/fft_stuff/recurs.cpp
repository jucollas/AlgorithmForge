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

constexpr int countr_zero_constexpr(unsigned int n) { int x = 0; while (!(n & (1 << x))) x++; return x; }
struct fft_info{
	static constexpr mint root =3;
	static constexpr int rank2=countr_zero_constexpr(mint::mod-1);
	vector<mint> wsq,iwsq;
	fft_info(){
		mint w=root.pow((mint::mod-1)>>rank2);
		wsq.resize(rank2+1);rep(i,rank2,-1){
			wsq[i]=w;w*=w;
		}
		w=wsq.back().inverse();
		iwsq.resize(rank2+1);rep(i,rank2,-1){
			iwsq[i]=w;w*=w;
		}
		
	}
};

#define INVERSE_POWER
void recursive_fft( fps &up,bool invert ){
	const int n=sz(up);	
	if(n<=1)return;
	
	fps odd(n/2),even(n/2);
	rep(i,0,n/2){
		odd[i] =up[i*2+1];
		even[i]=up[i*2];
	}
	recursive_fft(odd,invert);recursive_fft(even,invert);
	
	
	static const fft_info info;
#ifndef INVERSE_POWER
	const mint wn=(invert)?info.iwsq[ilog2(n)]:info.wsq[ilog2(n)];
#else
	const mint wn=info.wsq[ilog2(n)];
#endif
	mint w=1; rep(i,0,n/2){
		up[i]    =even[i]+w*odd[i];
		up[i+n/2]=even[i]-w*odd[i];
		w*=wn;
	}
}

void fft(fps&a,bool invert){
	/*fps brute(sz(a),0);
	static const fft_info info;
	mint wn=info.wsq[ilog2(sz(a))];
	if(invert)wn=wn.inverse();
	// assert (wn.pow(sz(a)/2)==mint(-1)); assert(wn.pow(sz(a))==mint(1));
	rep(i,0,sz(a)){
		mint wj=wn.pow(i+1),w=1;
		rep(j,0,sz(a)){
			brute[i]+=w*a[j];
			w*=wj;
		}
	}*/
	// idebug(a);
	recursive_fft(a,invert);
	// idebug(a);
	
	///this allows me to work identically to the other ffts I have to compare
	// reverse(all(a));reverse(a.begin()+1,a.end());
	
	// idebug(a);idebug(brute);
	
	if(invert){
#ifndef INVERSE_POWER
#else
		reverse(a.begin()+1,a.end());
#endif

        mint n_1 = mint(sz(a)).inverse();
        for (mint & x : a)x*=n_1;
		// idebug(a);
	}else{
		return;
		fps tmp=a;recursive_fft(tmp,1);
		idebug(tmp);
		mint n_1 = mint(sz(a)).inverse();
		// debug(n_1,sz(a),n_1*mint(sz(a)));
        for (mint & x : tmp)x*=n_1;
		idebug(tmp);
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
	// debug(nm);
	fft(A,0);fft(B,0);
	rep(i,0,sz(A))A[i]*=B[i];
	fft(A,1);
	// idebug(A);
	p_trunc(A,nm);
	// idebug(A);
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
	else rep(i,0,n)scanf("%d",&A[i]);
	
	fps B(n);
	if(ios)rep(i,0,n)cin>>B[i];
	else rep(i,0,n)scanf("%d",&B[i]);
	
	p_trunc(A,n); p_trunc(B,n);
	
	auto start=take_time();
	p_mult(B,A);
	auto milisec=get_durat(start);
	cerr << "recuss time " << milisec << '\n';
	if(ios)rep(i,0,2*n)cout << (i<sz(B)?B[i].vl:0) << ' ';
	else rep(i,0,2*n)printf("%d ",(i<sz(B)?B[i].vl:0));
	cout << '\n';
	
	return 0;
}