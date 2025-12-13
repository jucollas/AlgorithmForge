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

// mt19937_64 rng_64( chrono::steady_clock::now().time_since_epoch().count() );
mt19937_64 rng_64( 6546138465 );
constexpr int ilog2( int num ) { return 8*sizeof(int) - __builtin_clz( num ) - 1; }
template<typename tpow> constexpr tpow mpow(tpow x,lint e,tpow m){tpow res=1;while(e){if(e&1ll)res=(res*1ll*x)%m;e>>=1;x=(x*1ll*x)%m;}return res;}

const int mod = 998244353;

/*
Author: Oscar Vargas Pabon

It is though to work modulo primes, so inv (.inv,/,/=) and pow may not
	work properly otherwise.

I assume from my template :
inv :: int mpow(int x,int e,int m){int res=1;while(e){if(e&1)res=(res*1ll*x)%m;e>>=1;x=(x*1ll*x)%m;}return res;}	

Tested in testing/test_alghelp.cpp and in fft/ntt stuff
*/

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

const int template_limit = 1e6;
int a[template_limit], b[template_limit];

void primitive_root(int mdd){
	//tests for g such that $\forall_{p|(md-1)} g^{(md-1)/p}=1(mod md)$
	vector<int> dec;int md=mdd-1;
	int i=2;while(i*i<=md){
		if(md%i==0){
			dec.pb(i);
			while(md%i==0)md/=i;
		}
		++i;
	} if(md>1)dec.pb(md);
	auto pcond=[&](mint g)->bool{
		bool rs=1;for(int p:dec)rs&=!(g.pow((mdd-1)/p)==1); return rs;
	};
	i=2;while(!pcond(i))++i;
	std::cerr << i << " _ primitive root" << endl;
}


struct fft_info{
	static const int g=5;
	vector<int> dec;
	vector<mint> root,iroot;
	fft_info(){
		// factorization
		int md=mod-1; for(int ind=2;ind*ind<=md;++ind)while(md%ind==0){
			md/=ind;dec.push_back(ind);
		}if(md>1)dec.push_back(md);
		// shuffle(all(dec),rng_64);// just to really test it
		sort(all(dec));// for general cases???
		// idebug(dec);
		root.resize(dec.size()+1); iroot.resize(dec.size()+1);
		root.back()=g; iroot.back()=mint(g).inv();
		rep(i,sz(dec),0){
			root[i-1]=root[i].pow(dec[i-1]);
			iroot[i-1]=iroot[i].pow(dec[i-1]);
			// assert(root[i]*iroot[i]==1);
		}
	}
};

vector<mint> fft(const vector<mint> &A,int ind,bool invert){
	if(ind<0)return A;
	static fft_info info;
	
	const int k=info.dec[ind],n=sz(A);
	// debug(k,n,ind);
	vector<vector<mint>> ys(k);rep(i,0,k){
		vector<mint> as(n/k);rep(j,0,n/k){
			// debug(n/k,i,j,j*(n/k)+i);
			as[j]=A[j*k+i];
		}
		// debug(i);
		ys[i]=fft(as,ind-1,invert);
	}
	// debug("=???");
	vector<mint> y(n,0);
	// mint wr=(invert)?info.iroot[ind]:info.root[ind],wn,wf,w;wn=w.pow(n/k);
	// w -> wf^i ; wf -> w^{j*(n/k)+i} ; wn -> w_n ; wnj -> wn^j ; wnk -> wn^{n/k}
	mint w,wf,wn=(invert)?info.iroot[ind]:info.root[ind],wnj=1,wnk=wn.pow(n/k),wnki;
	// no, my method has a flawwww;  ; ; ; it works in 
	// T(n)= k T(n/k) + O(n)*O(k)
	// debug(wn,wn.pow(n));
	rep(i,0,k){
		wnki=1;
		rep(j,0,n/k){
			w=1; wf=wnj*wnki;
			// debug(w,wf,i,j);
			rep(l,0,k){
				y[j*k+i]+=w*ys[l][j];
				w*=wf;
			}
			wnki*=wnk;
		}
		wnj*=wn;
	}
	// idebug(y);
	return y;
}

vector<mint> naive_fft(vector<mint> A,int ind,bool invert){
	static fft_info info;
	mint wn=(invert)?info.iroot[ind]:info.root[ind];
	const int n=sz(A);
	vector<mint> C(n,0);
	rep(i,0,n)rep(j,0,n)C[i]+=wn.pow(i*j)*A[j];
	return C;
}

vector<mint> convolution(vector<mint> A,vector<mint> B){
	static fft_info info;
	const int mx=sz(A)+sz(B);
	
	int ind=-1,vl=1;while(vl<mx){
		// debug(mx,vl);
		vl*=info.dec[++ind];
	}
	// debug(mx,vl,ind);idebug(info.dec);
	A.resize(vl,0);B.resize(vl,0);
	
	vector<mint> aa=A,bb=B;
	A=fft(A,ind,0);
	// B=A;
	// mint w=info.root[ind];
	// vector<mint> C(vl,0);rep(i,0,vl)rep(j,0,vl){
		// C[i]+=w.pow(i*j)*B[j];
	// }
	// mint tmp=0;rep(i,0,vl)tmp+=B[i];debug(tmp);
	// idebug(A);idebug(C);idebug(B);
	
	// return {-1};
	B=fft(B,ind,0);
	aa=naive_fft(aa,ind,0);
	bb=naive_fft(bb,ind,0);
	idebug(A);idebug(B);idebug(aa);idebug(bb);
	
	rep(i,0,vl)A[i]*=B[i];
	aa=naive_fft(A,ind,1);
	A=fft(A,ind,1);
	
	idebug(A);idebug(aa);
	
	for(mint &ac:A)ac/=mint(vl);
	return A;
}
vector<mint> naive_convo(vector<mint>A,vector<mint> B){
	vector<mint>C(sz(A)+sz(B),0);rep(i,0,sz(A))rep(j,0,sz(B)){
		C[i+j]+=A[i]*B[j];
	}
	return C;
}

void solve() {
	int n;cin>>n;
	vector<mint> A(n),B(n);
	rep(i,0,n)cin>>A[i];
	rep(i,0,n)cin>>B[i];
	// I think im having some issue with the roots I choose
	vector<mint> C1=convolution(A,B),C2=naive_convo(A,B);
	if(C1!=C2){
		idebug(A);idebug(B);
		debug(A.front()*B.front());
		idebug(C1);idebug(C2);
	}else debug("ok");
}

int32_t main(){
	ios_base::sync_with_stdio(false);
    cin.tie(NULL);
	cout << setprecision(12) << fixed;
	
	// primitive_root(mod);
	
    int t = 2;
    // cin >> t; ++t;
    while ( --t ) {
		solve();
    }
	return 0;
}

