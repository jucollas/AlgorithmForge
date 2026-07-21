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
    // #define cerr for(;false;) cerr
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

const int mod=998244353;
#include"../../../modulo_int.cpp"
#include"../../../formalPowerSeries.cpp"
template<typename tint>vector<tint>vgen(int n){
	vector<tint> res(n);for(tint&ac:res)ac=lint(rng_64());
	return res;
}

template<typename tfps> std::vector<tfps> pow_proj(int k,
		const std::vector<tfps>&F,const std::vector<tfps>&W){
	assert(F.size()==W.size());
	if(F.empty()||k==0)return std::vector<tfps>(k,tfps(0));
	const int lgi=ilog2(W.size()-1)+1,nm=1<<lgi;
	std::vector<tfps> G(2*nm),Q(2*nm),GG(4*nm),QQ=GG;
	rep(i,0,nm)G[nm-i-1]=i<int(W.size())?W[i]:tfps(0),Q[i]=i<int(F.size())?-F[i]:tfps(0);
	const tfps c=F[0];Q[0]=0; // https://noshi91.hatenablog.com/entry/2023/12/10/163348
	std::vector<tfps>wrot(2*nm),wrotpw=wrot;{
		const tfps wlen=internal::rank_fft_root<tfps>[lgi+2];
        wrot[0]=1;for(int i=1;i<2*nm;++i)wrot[i]=wrot[i-1]*wlen;
        for(int i=0;i<2*nm;++i){ // bit-reversal
            int x=0;for(int e=0;e<lgi+1;++e)x|=(i>>e&1)<<(lgi-e);
            if(x<i)swap(wrot[i],wrot[x]);
        } for(int i=0;i<2*nm;++i)wrotpw[i]=wrot[i].pow(2*nm),wrot[i]=(wrot[i]+wrot[i]).inv();
    }int n=nm,m=1; while(n>1){
		for(int i=0;i<2*nm;++i)GG[i]=G[i],QQ[i]=Q[i];
		for(int i=2*nm;i<4*nm;++i)GG[i]=QQ[i]=0;

		internal::fft(GG,0); internal::fft(QQ,0);
		for(int i=0;i<2*nm;++i)G[i]=(
			(GG[i<<1]*QQ[i<<1|1]-GG[i<<1|1]*QQ[i<<1]+(GG[i<<1]-GG[i<<1|1])*wrotpw[i])
			)*wrot[i];
		for(int i=0;i<2*nm;++i)Q[i]=(QQ[i<<1]*QQ[i<<1|1])+(QQ[i<<1]+QQ[i<<1|1])*wrotpw[i];
		internal::fft(Q,1); internal::fft(G,1);
		for(int i=0;i<2*m;++i)for(int j=n/2;j<n;++j)G[i*n+j]=Q[i*n+j]=0;
		n>>=1;m<<=1;
	} for(int i=0;i<m;++i)G[i]=G[i*2];
	for(int i=0;i<m/2;++i)swap(G[i],G[m-i-1]);
	for(int i=m;i<2*m;++i)G[i]=0; if(c!=tfps(0)){ 
		std::vector<tfps> fc(k),ifc(k);//factorial/inverse-factorial
		fc[0]=1; for(int i=1;i<k;++i)fc[i]=fc[i-1]*tfps(i);
		ifc[k-1]=fc[k-1].inv();for(int i=k-1;i>0;--i)ifc[i-1]=ifc[i]*tfps(i);
		const int h=1<<(ilog2(2*k-1)+1); Q.resize(h);G.resize(h);
		Q[0]=1;for(int i=1;i<k;++i)Q[i]=Q[i-1]*c;
		for(int i=0;i<k;++i)Q[i]*=ifc[i],G[i]*=ifc[i];
		for(int i=k;i<h;++i)G[i]=Q[i]=0;
		internal::fft<tfps>(Q,0);internal::fft<tfps>(G,0);
		for(int i=0;i<h;++i)G[i]*=Q[i];
		internal::fft<tfps>(G,1); G.resize(k);
		for(int i=0;i<k;++i)G[i]*=fc[i];	
	} else G.resize(k,tfps(0));
	return G; }
vector<mint> brute(int m,const fps&F,const vector<mint>&r_W){
	const int n=r_W.size();fps W=r_W;reverse(all(W.F));
	fps B;rep(i,0,m)B[i]=(W*F.pow(i,n))[n-1];
	B.F.resize(m);
	return B.F;
}
template<typename tfps>
FormalPowerSeries<tfps> compositional_inverse(const FormalPowerSeries<tfps>&F,int n){
	// n|x^n|H(F(x))=|x^{-1}|H'(x)*G^{-n}(x) ; where G(F(x))=F(G(x))=x; H(x)=x^i for this
	assert(F[0]==tfps(0)&&F[1]!=tfps(0));//n|x^n|F(x)^i=|x^{n-i}|i(x/G(x))^n
	std::vector<tfps> W(n,tfps(0)),FF=F.F;W[n-1]=1;FF.resize(n,tfps(0));
	const tfps c=FF[1],ic=c.inv(); // to ensure there is a nth root
	for(int i=0;i<n;++i)FF[i]*=ic;
	std::vector<tfps> pp=pow_proj(n,FF,W);
	fps G(n-1);for(int i=1;i<n;++i)G[n-1-i]= pp[i]*tfps(n-1)*tfps(i).inv();
	G=G.pow(-tfps(n-1).inv(),n-1)<<1;
	tfps cc=1;for(int i=0;i<n;++i)G[i]*=cc,cc*=ic;
	return G;
}
template<typename tfps>
FormalPowerSeries<tfps> brute_composition(const FormalPowerSeries<tfps>&F,const FormalPowerSeries<tfps>&G,int n){
	//H=F(G(x))
	fps H;rep(i,0,n)H+=G.pow(i,n+1)*F[i];
	return H.trunc(n);
}
void solve() {
	int type=0;
	// type=1;
	if (type==1){
		const int n=1<<7,m=n,iter=100;
		rep(tm,0,iter){
			fps F=vgen<mint>(n);F[0]=0;
			// F={0, 48, 899930270 ,355122750 ,885751026, 105435380 ,936068561 ,694651622};
			fps G=compositional_inverse(F,m);
			fps FG=brute_composition(F,G,m);
			fps obj={0,1};
			if(FG!=obj){
				debug("bad_stuff");
				idebug(F);idebug(G);idebug(FG);
				assert(0);
			}
			// exit(0);
		}debug("Lo_logramos");
	} else {
		const int n=1<<17,m=n,iter=5;
				rep(tm,0,iter){
			fps F=vgen<mint>(n);F[0]=0;
			fps G=compositional_inverse(F,m);
		}debug("Lo_logramos");
	}
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