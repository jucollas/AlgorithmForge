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

mt19937 rng( chrono::steady_clock::now().time_since_epoch().count() );
constexpr int ilog2( int num ) { return 8*sizeof(int) - __builtin_clz( num ) - 1; }
constexpr int mpow(int x,int e,int m){int res=1;while(e){if(e&1)res=(res*1ll*x)%m;e>>=1;x=(x*1ll*x)%m;}return res;}

const int mod= 998244353;
/*
Autor:Oscar Vargas Pabon

mulo Integer
makes stuff a bit easier

pow assumes mpow from my template
*/

/*
Autor:Oscar Vargas Pabon

modulo Integer, makes stuff a bit easier

Note that pow and inv may not work well if the modulo is not prime

pow assumes mpow from my template
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
};

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
	constexpr montgomery_int( int v)noexcept:vl(transform(v>=0?(v<int(m)?v:v%m):(m-v<m?m-v:(v%m)+m))){};
	constexpr montgomery_int(lint v)noexcept:vl(transform(v>=0?(v<int(m)?v:v%m):(m-v<m?m-v:(v%m)+m))){};
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
	bool operator ==(const  uint &ot)const{return (*this)==montgomery_int(ot);}
	bool operator ==(const  int &ot )const{return (*this)==montgomery_int(ot);}
	bool operator !=(const montgomery_int &ot)const{return vl!=ot.vl;}
	
	operator bool() const { return itransform(vl); }
	operator  int() const { return itransform(vl); }
	
	friend ostream &operator<<(ostream &os,const montgomery_int &ac){return os << itransform(ac.vl);}
	friend istream &operator>>(istream&is,montgomery_int &ac){int v;is>>v;ac=montgomery_int(v);return is;}	
};


// typedef montgomery_int<mod> mint;
typedef modulo_int<mod> mint;

int main(){
	const int n=1e6;
	rep(i,0,n){
		int u=rng()%mod,v=rng()%mod;
		if(u<0)u+=mod;
		if(v<0)v+=mod;
		
		mint mu=u,mv=v,tmp;
		
		mint msum=mu+mv; int sum=(u+v)%mod;
		if(sum<0)sum+=mod;
		tmp=mu;tmp+=mv;
		if(!(msum ==sum) || tmp!=msum )debug(msum,sum,u,v);
		
		mint mres=mu-mv; int res=((u-v)%mod+mod)%mod;
		if(res<0)res+=mod;
		tmp=mu;tmp-=mv;
		if(!(mres ==res) || tmp!=mres )debug(mres,res,u,v);
		
		mint mneg=-mu; int neg=mod-u;
		if(!(mneg ==neg)  )debug(mneg,neg,u,v);
		
		
		mint mmul=mu*mv;int mul=(u*1ll*v)%mod;
		if(mul<0)mul+=mod;
		tmp=mu;tmp*=mv;
		if(!(mmul ==mul) || tmp!=mmul )debug(mmul,mul,u,v);
		
		mint mypow=mu.pow(v);int pow=mpow(u,v,mod);
		if(pow<0)pow+=mod;
		if(!(mypow==pow))debug(mypow,pow,u,v);
		
		int uu=mu;if(uu!=u)debug(u,uu,mu);
		
		if(v==0)continue;
		mint mdiv=mu/mv;int div=(u*1ll*mpow(v,mod-2,mod))%mod;
		if(div<0)div+=mod;
		tmp=mu;tmp/=mv;
		if(!(mdiv ==div) || tmp!=mdiv )debug(mdiv,div,u,v);
	}
	return 0;
}