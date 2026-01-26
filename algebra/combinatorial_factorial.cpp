/*
Author:Oscar Vargas Pabon

This is currently untested

I have 3 versions:
	1. O(n+lgn) precomputing ; O(1) answer | TESTED
	2. O(n^2lgn) using DP
	3. O(k) using the (n!/k!)*(1/k!) trick
*/

template<typename tcmb>
struct Combi{
	vector<tcmb> fact,ifact;
	constexpr Combi(int n){
		fact.resize(n);ifact.resize(n);
		fact[0]=1;rep(i,1,n)fact[i]=fact[i-1]*tcmb(i);
		ifact[n-1]=fact[n-1].inv();
		rep(i,n-1,0)ifact[i-1]=ifact[i]*tcmb(i);
	}
	tcmb comb(int n,int k)const{
		// n choose k
		if(n<0||k<0||n<k)return 0;
		return fact[n]*ifact[k]*ifact[n-k];
	}
};const Combi<mint> cmb(1e5);

template<typename tcmb>
struct Combi{
	map<pair<int,int>,tcmb>mem;
	Combi()=default;
	tcmb comb(int n,int k){
		tcmb res; if(n<0||k<0||n<k)res=0;
		else if(n==0||k==n)res=1;
		else if(mem.count({n,k}))res=mem[{n,k}];
		else{
			res=comb(n-1,k)+comb(n-1,k-1);
			mem[{n,k}]=res;
		}
		return res;
	}
};Combi<mint> cmb;

template<typename tcmb>
tcmb comb(int n,int k){
	if(n<0||k<0||n<k)return 0;
	if(n-k<k)return comb<tcmb>(n,n-k);
	tcmb nk=1;rep(i,n-k+1,n+1)nk*=i;
	tcmb kf=1;rep(i,1,k+1)kf*=i;
	return nk/kf;
}