/* Author: Oscar Vargas Pabon

A bunch of scattered implementations regarding fps and ntt stuff.
	circular_convolution: $(F*G)_k=\sum_{i+j=k(mod n)}F[i]*G[j]$
	chirpz: $F(x^0),F(x^1),F(x^2),...,F(x^n)$
	inner_product: $(F . G)_k=\sum_i F[i]*G[i+k]$
	prefab_hand: exemplifies stuff like $F=F*G+1$
	fps sqrt: using Tonelli shanks to get a sqrt
*/
vector<mint> circular_convolution(vector<mint>F,vector<mint>G){
    // Inspiration from the toeplitz matrix part on 
    // "Computational frameworks for the fast fourier transform"
    const int n=F.size(),lgi=ilog2(n-1)+2; if( (n&(n-1))==0 ){
        internal::fft(F,0);internal::fft(G,0);
        rep(i,0,n)F[i]*=G[i];
        internal::fft(F,1);return F;
    } F.resize(1<<lgi,0);G.resize(1<<lgi,0);
    rep(i,0,n)F[(1<<lgi)-n+i]=F[i];
    internal::fft(F,0);internal::fft(G,0);
    rep(i,0,1<<lgi)F[i]*=G[i];
    internal::fft(F,1);
    F.resize(n); return F;
} vector<mint> chirpz(vector<mint>A,const mint x=internal::primitive_root(mod)){
    const int n=A.size(),lgi=ilog2(n-1)+1;
    vector<mint>B(2<<lgi,0),ixs(n);A.resize(2<<lgi,0);
    const mint ix=x.inv(); mint ixa=1,xa=1;
    ixs[0]=1; rep(i,1,n)ixs[i]=ixs[i-1]*ixa,ixa*=ix,A[i]*=ixs[i];
    B[0]=1; rep(i,1,2*n)B[i]=B[i-1]*xa,xa*=x;
    internal::transposed_fft(B,1);internal::fft(A,0);
    rep(i,0,2<<lgi)A[i]*=B[i]; // this shit is symmetrical, lol
    internal::transposed_fft(A,0);
    rep(i,0,n)A[i]*=ixs[i];
    A.resize(n); return A;
} fps inner_product(fps A,fps B){
	const int n=A.size(),lgi=ilog2(n-1)+2;
	A.F.resize(1<<lgi,0);B.F.resize(1<<lgi,0);
	internal::fft(B.F,0); internal::transposed_fft(A.F,1);
	rep(i,0,1<<lgi)A[i]*=B[i];
	internal::transposed_fft(A.F,0); return A.trunc(n);
} fps redu_inner_product( fps A, fps B ){
	// C_k=\sum_{k=i-j}A_iB_j
	const int n=sz(A); reverse(all(B));
	return ((A*B)>>(n-1)).trunc(n);
} fps prefab_hand(const fps &d,const int n){
	//This is taken from fps_24/N . There is actually a simpler solution in O(nlgn)
	// to prefab using the fact that exp(log(H)) leaves exp( \sum_{r>0}d_r\sum_{i>0}x^{ir}/i )
	// I include this sol to exemplify online-fft

	// Generates the first n terms of $\prod_{r>0} (1-x^r)^{-d_r}$
	fps D(n,0); // I assume d[0]=0; D[0]=0;
	rep(i,1,d.size())if(i<=n)rep(j,1,n/i+1) D[j*i]+=d[i]*mint(i);
	
	fps H(n+1,0); function<void(int,int)>cdq=[&](int l,int r){
		if(l>=r) H[l]=(l)?H[l]/mint(l):mint(1);
		else{ const int m=(l+r)/2; cdq(l,m);
			// Left-Hand H[l..m]
			fps lh(m-l+1);rep(i,l,m+1)lh[i-l]=H[i];
			
			// Small-D : only relevant D terms
			int szd=min<int>(d.size(),r-l+1);
			fps sd(szd);rep(i,0,szd)sd[i]=D[i];
			
			lh*=sd; // Contribution of lh to H(m..r]
			rep(i,m+1,r+1)if(sz(lh)>i-l)H[i]+=lh[i-l];
			
			cdq(m+1,r);
		}
	};cdq(0,n); return H;
} namespace internal{
template<typename tfps> int Tonelli_Shanks(tfps a) {
	// usado por sqrt cuando trabajo con enteros modulo algo
	//plagiado epicamente de https://judge.yosupo.jp/submission/270105
	const int mod=tfps::mod();
	if (int(a) < 2) return a;
	if (mpow<int>(a, (mod - 1) / 2, mod) != 1) return -1;
	if (mod % 4 == 3) return mpow<int>(a, (mod + 1) / 4, mod);
 
	tfps b = 3; if (mod != 998244353) {
		while (mpow<int>(b, (mod - 1) / 2, mod) == 1) {
			b=tfps(int(rng_64()%(mod-3)) + 2);
		}
	} int q = mod - 1,Q = 0;
	while ( !(q&1) ) Q++, q /= 2;
 
	tfps x = mpow<int>(a, (q + 1) / 2, mod);
	b = mpow<int>(b, q, mod);
 
	int shift = 2; while ( x*x != a) {
		tfps error= tfps(mpow<int>(a,mod-2,mod))*x*x;
		if (mpow<int>(error, 1 << (Q - shift), mod) != 1) x *= b;
		b *=b; ++shift;
	} return x; }
} template<typename tfps>
FormalPowerSeries<tfps> fps_sqrt(const FormalPowerSeries<tfps>&F,int n){
	// this impl only works on numbers modulo a prime
	// first n terms of G^2=F O(nlgn)
	const tfps i2 = tfps(2).inv();//mpow(2,mod-2,mod);
	
	int xi=0;while(xi<F.size()&&F[xi]==0)++xi;
	if(xi>=F.size())return {0};
	
	if(xi&1)return {};
	int vl=internal::Tonelli_Shanks<tfps>(F[xi]);if(vl==-1)return {};
	FormalPowerSeries<tfps> G={vl},ac;//sqrt{F[0]}
	//I assume Tonelli_Shanks returns -1 when it's impossible
	
	FormalPowerSeries<tfps> H=(*this)>>xi;H.trunc(n-xi);
	
	for(int e=2;e<(n-xi)*2;e<<=1){
		ac=G.inv(e)*H;ac.trunc(e);
		G+=ac; G*=i2; G.trunc(e);
	} return (G<<(xi/2)).trunc(n);
}