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
}template<typename tfps,long long N> inline constexpr tfps fft_iroot(){return fft_root<tfps,N>().inv();}
// template<class tfps> constexpr auto prec_rank_fft_iroot() {
// 	constexpr int rank=countr_zero_constexpr(tfps::mod()-1);
//     std::array<tfps,rank+1> iroot; // precompute the roots
//     iroot[rank]=rank_fft_root<tfps>[rank].inv();
//     for(int i=rank-1;i>=0;--i)iroot[i]=iroot[i+1]*iroot[i+1];
// 	return iroot;
// }template<class tfps> inline constexpr auto rank_fft_iroot = prec_rank_fft_iroot<tfps>();
// template<typename tfps,long long N> inline constexpr tfps fft_iroot(){
// 	if constexpr(N<=0ll||countr_zero_constexpr(N)>=int(rank_fft_iroot<tfps>.size()))return tfps(0);
// 	return rank_fft_iroot<tfps>[countr_zero_constexpr(N)];
// }