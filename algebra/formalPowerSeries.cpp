/*
Author: Oscar Vargas Pabon

Based on code by MarcosK and other people
	https://codeforces.com/contest/438/submission/340901913
and multiple blogs all around the place
	https://cp-algorithms.com/algebra/polynomial.html#inverse-series_1
	https://codeforces.com/blog/entry/56422
	https://codeforces.com/blog/entry/12513 - problem E

Tested in
	https://judge.yosupo.jp/problem/inv_of_formal_power_series
	https://judge.yosupo.jp/problem/exp_of_formal_power_series
	https://judge.yosupo.jp/problem/log_of_formal_power_series
	https://judge.yosupo.jp/problem/pow_of_formal_power_series
	https://judge.yosupo.jp/problem/sqrt_of_formal_power_series

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

    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;

        if (i < j)
            swap(a[i], a[j]);
    }

    for (int len = 2; len <= n; len <<= 1) {
        int wlen = invert ? root_1 : root;
        for (int i = len; i < root_pw; i <<= 1)
            wlen = (int)(1LL * wlen * wlen % mod);

        for (int i = 0; i < n; i += len) {
            int w = 1;
            for (int j = 0; j < len / 2; j++) {
                int u = a[i+j], v = (int)(1LL * a[i+j+len/2] * w % mod);
                a[i+j] = u + v < mod ? u + v : u + v - mod;
                a[i+j+len/2] = u - v >= 0 ? u - v : u - v + mod;
                w = (int)(1LL * w * wlen % mod);
            }
        }
    }

    if (invert) {
        int n_1 = inverse(n, mod);
        for (int & x : a)
            x = (int)(1LL * x * n_1 % mod);
    }
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
void p_square(fps &A){
	// A'=A*A in O(nlgn)
	int nm=A.size()*2;
	int lgi=ilog2(nm-1)+1;A.resize(1<<lgi);
	fft(A,0);
	for(int &ac:A)ac=(ac*1ll*ac)%mod;
	fft(A,1);
	p_trunc(A,nm);
}

void p_add(fps &F, int vl, int xi=0){
	// F+ vl*x^{xi}
	if(sz(F)<=xi)F.resize(xi+1);
	F[xi]+=vl;
}
void p_scale(fps &F,int vl){
	// F*vl
	rep(i,0,sz(F))F[i]=(F[i]*1ll*vl)%mod;
}
void p_add(fps &A, const fps &B, int sgn=1 ){
	// A + B*sgn ; I assume sgn\in\{-1,1\}
	A.resize(max(sz(A),sz(B)));
	rep(i,0,min(sz(A),sz(B)))A[i]=(A[i]+B[i]*sgn)%mod;
}
fps p_deriv(const fps &F){
	fps G(max(1,sz(F)-1));// G=D(F)
	rep(i,1,sz(F))G[i-1]=(F[i]*1ll*i)%mod;
	return G;
}
fps p_inte(const fps &F){
	fps G(sz(F)+1);//D(G)=F
	rep(i,0,sz(F))G[i+1]=(F[i]*1ll*mpow(i+1,mod-2,mod))%mod;
	return G;
}

fps p_inv(const fps &F,int n){
	assert(F[0]);// G*F=1
	fps G={mpow(F[0],mod-2,mod)};
	for(int e=2;e<2*n;e<<=1){
		fps ac=G;// gives a ~/2 speedup
		p_mult(ac,{F.begin(),F.begin()+min(sz(F),e)});
		rep(i,0,sz(ac))ac[i]=(mod-ac[i])%mod;
		
		p_add(ac,2); p_mult(G,ac);
		p_trunc(G,e,0);
	}
	p_trunc(G,n);
	return G;
}
fps p_log(const fps &F,int n){
	//first n terms of G=ln(F)=inte(D(F)/F)
	assert(F[0]==1);
	fps G=p_inv(F,n),dF=p_deriv(F);
	p_mult(G,dF); p_trunc(G,n-1);
	return p_inte(G);
}

fps p_exp(const fps &F,int n){
	assert(!F[0]);// first n terms of G=exp(F) 
	fps G={1},ac;
	for(int e=2;e<n*2;e<<=1){
		ac=p_log(G,e); ac.resize(max(sz(ac),min(e,sz(F))));
		rep(i,0,sz(ac))ac[i]=(
			(i<sz(F)?F[i]:0)- ac[i])%mod;
		// for some reason working with negatives
		// destroys something (NTT???)
		rep(i,0,sz(ac))if(ac[i]<0)ac[i]+=mod;
		p_add(ac,1);
		
		p_mult(G,ac); p_trunc(G,e,0);
	}
	p_trunc(G,n);
	return G;
}
fps p_pow(const fps &F,lint e,int n){
	// first n terms of G=F^e O(nlgn)
	if(!e)return {1};
	int xi=0;while(xi<sz(F)&&F[xi]==0)++xi;
	if( xi>=sz(F) || xi > n/e )return {0};
	
	int alp=mpow(F[xi],e%(mod-1),mod),ainv=mpow(F[xi],mod-2,mod);
	fps H(F.begin()+xi,F.end());p_scale(H,ainv);
	p_trunc(H,n-xi);
	
	int rx=xi*e;
	H=p_log(H,n-rx);
	p_scale(H,e%mod);
	H=p_exp(H,n-rx);
	
	fps G(n,0); H.resize(n-rx);
	rep(i,rx,n)G[i]=(H[i-rx]*1ll*alp)%mod;
	return G;
}
fps p_binpow( fps F, lint e, int n){
	// first n terms of G=F^e O(nlgnlge)
	fps G={1};while(e){
		if(e&1ll)p_mult(G,F);
		e>>=1;p_square(F);
		p_trunc(G,n);p_trunc(F,n);
	}
	return G;
}

int Tonelli_Shanks(int a, int mod) {
	//plagiado epicamente de https://judge.yosupo.jp/submission/270105
	
	if (a < 2) return a;
	if (mpow(a, (mod - 1) / 2, mod) != 1) return -1;
	if (mod % 4 == 3) return mpow(a, (mod + 1) / 4, mod);

	int b = 3;
	if (mod != 998244353) {
		while (mpow(b, (mod - 1) / 2, mod) == 1) {
			b=rng_64()%(mod-3) + 2;
		}
	}

	int q = mod - 1,Q = 0;
	while ( !(q&1)) Q++, q /= 2;

	int x = mpow(a, (q + 1) / 2, mod);
	b = mpow(b, q, mod);

	int shift = 2;
	while ((x *1ll* x) % mod != a) {
		int error = (((mpow(a, mod - 2, mod) *1ll* x) % mod) *1ll* x) % mod;
		if (mpow(error, 1 << (Q - shift), mod) != 1) {
			x = (x *1ll* b) % mod;
		}
		b = (b *1ll* b) % mod;
		++shift;
	}
	return x;
}


fps p_sqrt( const fps &F,int n){
	// first n terms of G^2=F O(nlgn)
	const int i2 = mpow(2,mod-2,mod);
	
	int xi=0;while(xi<sz(F)&&F[xi]==0)++xi;
	if( xi>=sz(F) )return {0};
	
	fps G={Tonelli_Shanks(F[xi],mod)},ac;//sqrt{F[0]}
	if(G[0]==-1||(xi&1))return {};
	//I assume Tonelli_Shanks returns -1 when it's impossible
	// debug(xi,G.front());
	fps H(F.begin()+xi,F.end()); p_trunc(H,n-xi);
	
	
	
	for(int e=2;e<(n-xi)*2;e<<=1){
		ac=p_inv(G,e);
		p_mult(ac,H);p_trunc(ac,e,0);
		p_add(G,ac);p_scale(G,i2);
		p_trunc(G,e);
	}
	swap(H,G);
	G.resize(sz(G)+xi/2);
	rep(i,0,xi/2)G[i]=0;
	rep(i,xi/2,sz(G))G[i]=H[i-xi/2];
	return G;
}
