	
// tuning_fps.cpp says it doesnt work
// my conjecture is that it has to do with the weird reorderings we are evading
// by using the transposed algorithm
	FormalPowerSeries<tfps> pow2(lint e,int n)const{
		// first n terms of G=F^e O(n(lgn+lge))
		assert(e>=0); if(!e)return {1};
		std::vector<tfps> G={1},xi=F;
		int tsz=1;while(tsz<n)tsz*=2;
		G.resize(tsz,0);xi.resize(tsz,0);
		internal::fft(G,0);internal::fft(xi,0);
		while(e){
			debug(e);
			if(e&1ll)rep(i,0,tsz)G[i]*=xi[i];
			e>>=1;  rep(i,0,tsz)xi[i]*=xi[i];
		}
		internal::fft(G,1);
		return FormalPowerSeries<tfps>(G).trunc(n);
	}