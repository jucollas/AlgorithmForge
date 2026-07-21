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
// #include"../../../bostanMori.cpp"
template<typename tint>vector<tint>vgen(int n){
	vector<tint> res(n);for(tint&ac:res)ac=lint(rng_64());
	return res;
}
template<typename tfps> tfps bostanMori2(std::vector<tfps>P,std::vector<tfps>Q,__uint64_t k){
    // computes |x^k|P/Q in time O(d*lgd*lgk)
    const int d=Q.size(),d2=1<<(ilog2(d-1)+2);
    P.resize(d2,tfps(0));Q.resize(d2,tfps(0)); while(k){ 

        internal::fft(Q,0);internal::fft(P,0);
    	for(int i=0;i<d2;++i)P[i]*=Q[i^1];	
    	for(int i=0;i<d2;i+=2)Q[i]=Q[i|1]=Q[i]*Q[i|1];
        internal::fft(Q,1);internal::fft(P,1);
        
        for(int i=0;i<d;++i)Q[i]=Q[i*2],P[i]=P[i*2+(k&1ll)];
        k>>=1; for(int i=d;i<d2;++i)Q[i]=P[i]=0;
    } return P[0]/Q[0];
}

template<typename tfps> tfps old_bostanMori(std::vector<tfps>P,std::vector<tfps>Q,__uint64_t k){
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
template<typename tfps> tfps bostanMori3(std::vector<tfps>P,std::vector<tfps>Q,__uint64_t k){
    // computes |x^k|P/Q in time O(d*lgd*lgk)
    const int d=Q.size(),dlgi=ilog2(d-1)+2,d2=1<<dlgi;
    P.resize(d2,tfps(0));Q.resize(d2,tfps(0));
    internal::fft(Q,0);
    vector<tfps>Q2(d2/2); while(k){
    	internal::fft(P,0);
    	for(int i=0;i<d2;++i)P[i]*=Q[i^1];	
        
    	for(int i=0;i<d2/2;++i)Q2[i]=Q[i]=Q[i<<1]*Q[i<<1|1];
    	internal::fft(Q2,1); // My first fft-doubling
    	const tfps wlen=internal::rank_fft_root<tfps>[dlgi];
    	tfps w=1;for(int i=0;i<d2/2;++i){
        	Q2[i]*=w; w*=wlen; // noishi is the goat
        }internal::fft(Q2,0);
        for(int i=d2/2;i<d2;++i)Q[i]=Q2[i-d2/2];

        internal::fft(P,1);
        for(int i=0;i<d;++i)P[i]=P[i*2+(k&1ll)]; // Q[i]=Q[i*2]
        k>>=1; for(int i=d;i<d2;++i)P[i]=0; // Q[i]
    } internal::fft(Q,1);
    return P[0]/Q[0];
}
template<typename tfps> tfps bostanMori(std::vector<tfps>P,std::vector<tfps>Q,__uint64_t k){
    const int d=Q.size(),dlgi=ilog2(d-1)+2,d2=1<<dlgi,d22=d2/2;
    P.resize(d2,tfps(0));Q.resize(d2,tfps(0)); // computes |x^k|P/Q in time O(d*lgd*lgk)
    // https://noshi91.hatenablog.com/entry/2023/12/10/163348
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

 void solve() {
	const int type=0;
	// const int type=1;
	if constexpr(type==1){
		const int n=1<<10,iter=100;
		rep(tm,0,iter){
			vector<mint> F=vgen<mint>(n),W=vgen<mint>(n);
			lint k=abs(lint(rng_64()));
			// k=abs(lint(rng_64()))%7;
            k=1;

			mint Bm=bostanMori(F,W,k);
			mint Bb=old_bostanMori(F,W,k);
			if(Bm!=Bb){
				idebug(F);idebug(W);
				debug(k,Bb,Bm);
				assert(Bb==Bm);
			}
			// exit(0);
		}debug("Lo_logramos");
	} else {
		const int n=1<<18,m=n,iter=5;
		rep(tm,0,iter){
			debug(tm);
			vector<mint> F=vgen<mint>(n),W=vgen<mint>(n);
			const lint k=abs(lint(rng_64()));
			mint Bm=bostanMori(F,W,k);
		}
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