/* Author: Oscar Vargas Pabon

A bunch of scattered implementations regarding fps and ntt stuff.
	circular_convolution: $(F*G)_k=\sum_{i+j=k(mod n)}F[i]*G[j]$
	chirpz: $F(x^0),F(x^1),F(x^2),...,F(x^n)$
	inner_product: $(F . G)_k=\sum_i F[i]*G[i+k]$
	prefab_hand: exemplifies stuff like $F=F*G+1$
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
} vector<mint> chirpz(vector<mint>A,const mint x=internal::primitive_root(mod) ){
    // some stuff adapted from * https://nyaannyaan.github.io/library/ntt/chirp-z.hpp
    // and * https://codeforces.com/blog/entry/83532
    const int n=A.size(),lgi=ilog2(n-1)+1;
    vector<mint> B(2<<lgi,0),ixs(n);
    const mint ix=x.inv(); mint ixa=1,xa=1;
    ixs[0]=1; rep(i,1,n)ixs[i]=ixs[i-1]*ixa,ixa*=ix,A[i]*=ixs[i];
    B[0]=1; rep(i,1,2*n)B[i]=B[i-1]*xa,xa*=x;
    reverse(A.begin(),A.end());A.resize(2<<lgi,0);
    internal::fft(A,0);internal::fft(B,0);
    rep(i,0,2<<lgi)A[i]*=B[i];
    internal::fft(A,1);A=vector<mint>(A.begin()+n-1,A.begin()+2*n-1);
    rep(i,0,n)A[i]*=ixs[i];
    return A;
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
}
