// Old inverse

FormalPowerSeries<tfps> inv(int n) const {
		assert(F[0]);// G*F=1
		FormalPowerSeries<tfps> G={F[0].inv()};//mpow(F[0],mod-2,mod)
		for(int e=2;e<2*n;e<<=1){
			FormalPowerSeries<tfps> ac=G;
			ac*= FormalPowerSeries({F.begin(),F.begin()+min<int>(F.size(),e)}); // gives a ~/2 speedup
			for(tfps &act:ac)act=-act;
			ac+={2}; G*=ac;
			G.trunc(e);
		} return G.trunc(n);
	}
emplate<int N>void rec_inv(tfps*FF,tfps*G)const{
		if constexpr(N==1)G[0]=F[0].inv();
		if constexpr(N<=1)return;
		rec_inv<N/2>(FF,G); constexpr int N2=N*2;
		const int k=N<int(F.size())?N:F.size(); 
		for(int i=0;i<k;++i)FF[i]=F[i];
		for(int i=k;i<N2;++i)FF[i]=0;
		internal::butterfly_rec<tfps,bool(0),N2>(FF);
		internal::butterfly_rec<tfps,bool(0),N2>(G);
		for(int i=0;i<N2;++i)G[i]=G[i]+G[i]-G[i]*G[i]*FF[i];
		internal::ibutterfly_rec<tfps,bool(0),N2>(G);
		for(int i=0;i<N ;++i)G[i]*=tfps(N2).inv();
		for(int i=N;i<N2;++i)G[i] =0;
	} template<int N=std::min<long long>(1<<25,(tfps::mod()-1)&(1-tfps::mod()))>
	FormalPowerSeries<tfps> inv(int n) const {
		if constexpr(N<1) return FormalPowerSeries<tfps>();
		constexpr int N2=N/2; if(N2<n){ assert(F[0]);
			std::vector<tfps> FF(N*2),GF(N*2,tfps(0));
			rec_inv<N>(FF.data(),GF.data());
			return FormalPowerSeries<tfps>(GF).trunc(n);
		} else return inv<N2>(n);
	}
	
// precomputing both roots and inverse-roots
template<class tfps> constexpr auto prec_rank_fft_root() {
	constexpr int rank=countr_zero_constexpr(tfps::mod()-1);
    std::array<tfps,rank+1> root; // precompute the roots
    root[rank]=tfps(primitive_root(tfps::mod())).pow((tfps::mod()-1)>>rank);
    for(int i=rank-1;i>=0;--i)root[i]=root[i+1]*root[i+1];
	return root;
}template<class tfps> inline constexpr auto rank_fft_root = prec_rank_fft_root<tfps>();
template<typename tfps,long long N> inline constexpr tfps fft_root(){
	if constexpr(N<=0ll||countr_zero_constexpr(N)>=int(rank_fft_root<tfps>.size()))return tfps(0);
	return rank_fft_root<tfps>[countr_zero_constexpr(N)];
} template<class tfps> constexpr auto prec_rank_fft_iroot() {
	constexpr int rank=countr_zero_constexpr(tfps::mod()-1);
    std::array<tfps,rank+1> iroot; // precompute the roots
    iroot[rank]=rank_fft_root<tfps>[rank].inv();
    for(int i=rank-1;i>=0;--i)iroot[i]=iroot[i+1]*iroot[i+1];
	return iroot;
}template<class tfps> inline constexpr auto rank_fft_iroot = prec_rank_fft_iroot<tfps>();
template<typename tfps,long long N> inline constexpr tfps fft_iroot(){
	if constexpr(N<=0ll||countr_zero_constexpr(N)>=int(rank_fft_iroot<tfps>.size()))return tfps(0);
	return rank_fft_iroot<tfps>[countr_zero_constexpr(N)];
}

// old bostanMori
mint bostanMori( lint k, fps P, fps Q ){
	// computes |x^k|P/Q in time O(d*lgd*lgk)
	const int d=Q.size(); while(k){ 
		fps nQ=Q;rep(i,0,nQ.size())if(i&1)nQ[i]=-nQ[i];
		P*=nQ; Q*=nQ;
		
		rep(i,0,d)Q[i]=Q[i*2];
		Q.trunc(d);
	
		// P[d*2] to prevent issues of P[ind1]=P[ind2] when sz(P)>=ind2
		P[d*2]; rep(i,0,d)P[i]=P[i*2+(k&1ll)];
		P.trunc(d); k>>=1;
	} return P[0]/Q[0]; // the method below gives a rough ~1/4 speedup
} template<typename tfps> tfps bostanMori(std::vector<tfps>P,std::vector<tfps>Q,__uint64_t k){
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