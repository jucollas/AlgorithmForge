/* Author: Oscar Vargas Pabon
Impl guided by
	* https://maspypy.com/fps-composition-and-compositional-inverse-part-1-compositional-inverse-and-power-projection
	* https://maspypy.com/fps-composition-and-compositional-inverse-part-2
	* https://noshi91.hatenablog.com/entry/2023/12/10/163348
	* "Power Series Composition in Near-Linear Time"
Power projection:
	Given F,W,k, calculate $R_i=\sum_{j\geq0}W_j|x^j|F(x)^i$ for $i\in 0,1,...,k-1$
Compositional inverse:
	Given F,n, calculate first n terms of G such that $F(G(x))=G(F(x))=x$
Composition:
	Given F,G,n, calculate first n terms of $F(G(x))$
Composition is the transposed problem of power projection.
*/
template<typename tfps> std::vector<tfps> pow_proj(
		const std::vector<tfps>&F,const std::vector<tfps>&W,int k){
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
	return G;
} template<typename tfps>
FormalPowerSeries<tfps> compositional_inverse(const FormalPowerSeries<tfps>&F,int n){
	// n|x^n|H(F(x))=|x^{-1}|H'(x)*G^{-n}(x) ; where G(F(x))=F(G(x))=x; H(x)=x^i for this
	assert(F[0]==tfps(0)&&F[1]!=tfps(0));//n|x^n|F(x)^i=|x^{n-i}|i(x/G(x))^n
	std::vector<tfps> W(n,tfps(0)),FF=F.F;W[n-1]=1;FF.resize(n,tfps(0));
	const tfps c=FF[1],ic=c.inv(); // to ensure there is a nth root
	for(int i=0;i<n;++i)FF[i]*=ic;
	std::vector<tfps> pp=pow_proj(FF,W,n);
	fps G(n-1);for(int i=1;i<n;++i)G[n-1-i]= pp[i]*tfps(n-1)*tfps(i).inv();
	G=G.pow(-tfps(n-1).inv(),n-1)<<1;
	tfps cc=1;for(int i=0;i<n;++i)G[i]*=cc,cc*=ic;
	return G;
} template<typename tfps> std::vector<tfps> composition(
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
