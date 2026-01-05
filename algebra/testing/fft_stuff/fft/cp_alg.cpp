#include <bits/stdc++.h>

typedef long long lint;

using namespace std;

#define debug(args...) { string _s = #args; replace(_s.begin(), _s.end(), ',', ' '); stringstream _ss(_s); istream_iterator<string> _it(_ss); raw_debug(_it, args);}
void raw_debug(istream_iterator<string> it) {cerr<<endl;}
template<typename T, typename... Args>
void raw_debug(istream_iterator<string> it, T a, Args... args) { cerr <<"<"<< *it << "->" << a << "> "; raw_debug(++it, args...); }
#define idebug(v) {cout<<'['<<#v<<']';for(const auto &el:v)cout << ' ' << el; cout << endl;}
#define adebug(ar,n) {cout<<'['<<#ar<<']';for(int i=0;i<n;++i)cout << ' ' << ar[i]; cout << endl;}

#define rep(i,strt,end) for(int i = strt ; i !=int(end) ; (int(strt)<int(end))?++i:--i )
// #define rall(vec) vec.rbegin(), vec.rend()
// #define all(vec) vec.begin(), vec.end()
#define sz(vec) int(vec.size())
// #define pb push_back
#define pob pop_back
// #define pf push_front
// #define pof pop_front

mt19937_64 rng_64( chrono::steady_clock::now().time_since_epoch().count() );
int ilog2( int num ) { return 8*sizeof(int) - __builtin_clz( num ) - 1; }
int mpow(int x,int e,int m){int res=1;while(e){if(e&1)res=(res*1ll*x)%m;e>>=1;x=(x*1ll*x)%m;}return res;}


/*
Author: Oscar Vargas Pabon

Based on code by MarcosK
	https://codeforces.com/contest/438/submission/340901913
and multiple blogs all around the place
	like the cp-algo one
The first version of this impl is in atcoder fps_24 A

I assume from my template:
everything     :: #define rep(i,strt,end) for(int i = strt ; i !=int(end) ; (int(strt)<int(end))?++i:--i )
everything     :: #define sz(vec) int(vec.size())
p_trunc        :: #define pob pop_back
p_mult,p_square:: int ilog2( int num ) { return 8*sizeof(int) - __builtin_clz( num ) - 1; }
modInverse     :: int mpow(int x,int e,int m){int res=1;while(e){if(e&1)res=(res*1ll*x)%m;e>>=1;x=(x*1ll*x)%m;}return res;}	
TonelliShanks  :: mt19937_64 rng_64( chrono::steady_clock::now().time_since_epoch().count() );
*/

typedef int tfps;
typedef vector<tfps> fps;

/* START OF NTT */
const int mod = 998244353;
const int root = mpow(3,119,mod);
const int root_1 = mpow(root,mod-2,mod);
const int root_pw = 1 << 23;

int inverse(int vl,int md){return mpow(vl,md-2,md);}

void fft(fps & a, bool invert) {
    int n = a.size();
	// vector<int> nm(n);rep(i,0,n)nm[i]=i;
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;

        if (i < j)
            swap(a[i], a[j]);//,swap(nm[i],nm[j]);
    }
	// idebug(nm);
	// idebug(a);
    for (int len = 2; len <= n; len <<= 1) {
        int wlen = invert ? root_1 : root;
        for (int i = len; i < root_pw; i <<= 1)
            wlen = (int)(1LL * wlen * wlen % mod);
		// debug(len);
        for (int i = 0; i < n; i += len) {
            int w = 1;//int cnt=0;
            for (int j = 0; j < len / 2; j++) {
                int u = a[i+j], v = (int)(1LL * a[i+j+len/2] * w % mod);
				// debug(i+j,i+j+len/2,nm[i+j],nm[i+j+len/2], cnt,w,wlen);
				// debug((u+v)%mod,((u-v)%mod+mod)%mod,u,v);
                a[i+j] = u + v < mod ? u + v : u + v - mod;
                a[i+j+len/2] = u - v >= 0 ? u - v : u - v + mod;
                w = (int)(1LL * w * wlen % mod);// ++cnt;
            }
        }
    }

    if (invert) {
        int n_1 = inverse(n, mod);
        for (int & x : a)
            x = (int)(1LL * x * n_1 % mod);
    }
	// else idebug(a);
}
/* End of NTT */

void p_trunc(fps &F, int n,bool elim_0=1){
	F.resize(max(1,min(sz(F),n)));
	if(elim_0)while(sz(F)>1&&F.back()==0)F.pob();
}

void p_mult(fps &A, fps B){
	// A'=A*B in O(nlgn)
	int nm=A.size()+B.size();
	int lgi=ilog2(nm-1)+1;A.resize(1<<lgi);B.resize(1<<lgi);
	fft(A,0);fft(B,0);
	rep(i,0,sz(A))A[i]=(A[i]*1ll*B[i])%mod;
	fft(A,1);
	p_trunc(A,nm);
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
	cerr << "cpalg time " << milisec << '\n';
	if(ios)rep(i,0,2*n)cout << (i<sz(B)?int(B[i]):0) << ' ';
	else rep(i,0,2*n)printf("%d ",(i<sz(B)?int(B[i]):0));
	cout << '\n';
	
	return 0;
}

/*int main(){
	ios_base::sync_with_stdio(false);
    cin.tie(NULL);
	cout << setprecision(12) << fixed;
	int n;cin>>n;
	
	// lint pm=0;cin>>pm;
	
	fps A(n);rep(i,0,n)cin>>A[i];
	fps B(n);rep(i,0,n)cin>>B[i];
	p_trunc(A,n); p_trunc(B,n);
	
	auto start=take_time();
	p_mult(B,A);
	auto milisec=get_durat(start);
	cerr << "cp algo time " << milisec << '\n';
	
	
	// int i1=mpow(2,pm%mod,mod),i2=mpow(2,pm%(mod-1ll),mod);
	// debug(i1,i2);
	
	// fps B=p_inv(A,n);
	// fps B=p_log(A,n);
	// fps B=p_exp(A,m,n);
	//fps B=p_pow(A,pm,n);
	// fps C=p_binpow(A,pm,n); idebug(C);
	// fps D=p_binpow(A,pm%(mod-1),n); idebug(D);
	// fps E=p_binpow(A,mod+1,n);idebug(E);
	
	// fps B = p_sqrt(A,n);if(B.empty()){cout << "-1\n";return 0;}
	// fps C=B;p_square(C);p_trunc(C,n);idebug(C);
	rep(i,0,2*n)cout << (i<sz(B)?B[i]:0) << ' ';
	cout << '\n';
	
	return 0;
}*/