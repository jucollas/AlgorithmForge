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
// #pragma GCC optimize("Ofast")
// #define NDEBUG
#include <bits/stdc++.h>
#include <cassert>

typedef long long lint;
typedef __int128_t int128;
using namespace std;
#ifdef OSVARP
    #include<sys/resource.h>
#else
    #define cerr for(;false;) cerr
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

#include"clean_impl.cpp"

template<typename tfps> std::vector<tfps> composition(
        const std::vector<tfps>&H,const std::vector<tfps>&F,int k){
    assert(F.size()==H.size()&&F[0]==tfps(0));//first k terms of H(F(x))
    if(F.empty()||k==0)return std::vector<tfps>(k,tfps(0));
    const int lgi=ilog2(H.size()-1)+1,nm=1<<lgi;
    std::vector<tfps> G(2*nm),Q=G,GG(4*nm);
    for(int i=0;i<nm;++i)Q[i]=i<int(F.size())?-F[i]:tfps(0);
    std::vector<tfps>wrot(2*nm),wrotpw=wrot;{
        const tfps wlen=internal::rank_fft_root<tfps>[lgi+2];
        wrot[0]=1;for(int i=1;i<2*nm;++i)wrot[i]=wrot[i-1]*wlen;
        for(int i=0;i<2*nm;++i){ // bit-reversal
            int x=0;for(int e=0;e<lgi+1;++e)x|=(i>>e&1)<<(lgi-e);
            if(x<i)swap(wrot[i],wrot[x]);
        } for(int i=0;i<2*nm;++i)wrotpw[i]=wrot[i].pow(2*nm),wrot[i]=(wrot[i]+wrot[i]).inv();
    } auto calc=[&](auto rec,int n,int m)->void{
        if(n<=1)return;
        std::vector<tfps> QQ(4*nm,tfps(0));
        for(int i=0;i<2*nm;++i)QQ[i]=Q[i];
        internal::fft(QQ,0);
        for(int i=0;i<2*nm;++i)Q[i]=(QQ[i<<1]*QQ[i<<1|1])+(QQ[i<<1]+QQ[i<<1|1])*wrotpw[i];
        internal::fft(Q,1);
        for(int i=0;i<2*m;++i)for(int j=n/2;j<n;++j)Q[i*n+j]=0;

        rec(rec,n>>1,m<<1);

        for(int i=0;i<2*m;++i)for(int j=n/2;j<n;++j)G[i*n+j]=0;
        internal::transposed_fft(G,1);
        for(int i=0;i<2*nm;++i){
            GG[i<<1  ]= G[i]*(QQ[i<<1|1]+wrotpw[i])*wrot[i];
            GG[i<<1|1]=-G[i]*(QQ[i<<1  ]+wrotpw[i])*wrot[i];
        } internal::transposed_fft(GG,0);
        for(int i=0;i<2*nm;++i)G[i]=GG[i];
    };
    for(int i=nm-1;i>=0;--i)G[i*2]=nm-i-1<int(H.size())?H[nm-i-1]:tfps(0);
    calc(calc,nm,1);
    G.resize(nm);std::reverse(G.begin(),G.end());
    return G;
}


void solve() {
	int type=0;
    // type=1;
    if (type==1){
        const int n=1<<9,m=n,iter=100;
        rep(tm,0,iter){debug(tm);
            fps F=vgen<mint>(n),G=vgen<mint>(n);
            G[0]=0;
            vector<mint> FGm=composition(F.F,G.F,m);
            fps FGb=brute_composition(F,G,m);
            if(fps(FGm)!=FGb){
                debug("bad_stuff");
                idebug(F);idebug(G);idebug(FGm);idebug(FGb);
                assert(0);
            }
            // exit(0);
        }debug("Lo_logramos");
    } else {
        const int n=1<<17,m=n,iter=5; rep(tm,0,iter){debug(tm);
            vector<mint> F=vgen<mint>(n),G=vgen<mint>(n);G[0]=0;
            vector<mint> FGm=composition(F,G,m);
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