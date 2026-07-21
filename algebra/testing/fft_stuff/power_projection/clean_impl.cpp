const int mod=998244353;
#include"../../../modulo_int.cpp"
#include"../../../formalPowerSeries.cpp"
template<typename tint>vector<tint>vgen(int n){
	vector<tint> res(n);for(tint&ac:res)ac=lint(rng_64());
	return res;
}



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


vector<mint> brute_pow_proj(int m,const fps&F,const vector<mint>&r_W){
	const int n=r_W.size();fps W=r_W;reverse(all(W.F));
	fps B;rep(i,0,m)B[i]=(W*F.pow(i,n))[n-1];
	B.F.resize(m);
	return B.F;
}fps brute_composition(const fps&F,const fps&G,int n){
	//H=F(G(x))
	fps H;rep(i,0,n)H+=G.pow(i,n+1)*F[i];
	return H.trunc(n);
}