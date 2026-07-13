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
template<typename tpow,typename texp=lint> constexpr tpow mpow(tpow x,lint e,tpow m){tpow res=1;while(e){if(e&1ll)res=(res*texp(1)*x)%m;e>>=1;x=(x*texp(1)*x)%m;}return res;}

const int mod= 998244353;
#define MONTGOMERY
#ifdef MONTGOMERY
	#include"../montgomery_space.cpp"
#else
	#include"../modulo_int.cpp"
#endif


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