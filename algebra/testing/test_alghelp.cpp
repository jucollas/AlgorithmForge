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


/*
Autor:Oscar Vargas Pabon

Modulo Integer
makes stuff a bit easier
*/
struct modulo_int{
	static const int mod= 998244353; static_assert(mod>0);
	
	int vl;
	constexpr modulo_int()noexcept:vl(0){};
	constexpr modulo_int( int v)noexcept:vl(v>=0?v%mod:(v%mod)+mod){};
	constexpr modulo_int(lint v)noexcept:vl(v>=0?v%mod:(v%mod)+mod){};
	constexpr modulo_int(unsigned long long v)noexcept:vl(v%mod){};
	
	modulo_int &operator +=(const modulo_int &ot){ vl+=ot.vl; if(vl>=mod)vl-=mod; return *this; }
	modulo_int  operator + (const modulo_int &ot)const{ return modulo_int(*this)+=ot; }
	modulo_int &operator -=(const modulo_int &ot){ vl-=ot.vl; if(vl<0)vl+=mod; return *this; }
	modulo_int  operator - (const modulo_int &ot)const{ return modulo_int(*this)-=ot; }
	modulo_int &operator *=(const modulo_int &ot){ vl=(vl*1ll*ot.vl)%mod; return *this; }
	modulo_int  operator * (const modulo_int &ot)const{ return modulo_int(*this)*=ot; }
	modulo_int &operator /=(const modulo_int &ot){ (*this)*=ot.inverse(); return *this; }
	modulo_int  operator / (const modulo_int &ot)const{ return modulo_int(*this)/=ot; }
	
	modulo_int inverse()const{return modulo_int(mpow(vl,mod-2,mod));}//Fermats little theorem
	modulo_int operator -()const {return modulo_int(-vl);}
	modulo_int pow(lint e)const{return modulo_int(mpow(vl,e%(mod-1),mod));}
	
	bool operator ==(const modulo_int &ot)const{return vl==ot.vl;} 
	bool operator ==(const  int &ot)const{return vl==ot;   }
	bool operator !=(const modulo_int &ot)const{return vl!=ot.vl;}
	
	operator bool() const { return vl; }
	operator  int() const { return vl; }
	
	friend ostream &operator<<(ostream &os,const modulo_int &ac){return os << ac.vl;}
	friend istream &operator>>(istream&is,modulo_int &ac){int v;cin>>v;ac=modulo_int(v);return is;}	
}; //const int mod=modulo_int::mod; typedef modulo_int mint;

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
		// debug(x,"enestaver",q,m,res);
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
	
	montgomery_int operator -()const {return mod-vl;}
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
		if(v==0)continue;
		mint mdiv=mu/mv;int div=(u*1ll*mpow(v,mod-2,mod))%mod;
		if(div<0)div+=mod;
		tmp=mu;tmp/=mv;
		if(!(mdiv ==div) || tmp!=mdiv )debug(mdiv,div,u,v);
		mint mypow=mu.pow(v);int pow=mpow(u,v,mod);
		if(pow<0)pow+=mod;
		if(!(mypow==pow))debug(mypow,pow,u,v);
	}
	
	rep(i,0,n){
		int u=rng()%mod;if(u<0)u+=mod;
		mint mu=u;
		int uu=mu;
		assert(uu==u);
	}
	
	return 0;
}