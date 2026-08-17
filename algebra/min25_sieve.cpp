/* Author: Oscar Vargas Pabon
This is mostly to guide myself. I assume there exists some function f 
	that is the multiplicative function we are using. There is also some 
	weird g that must also be multiplicative and will be used to help in the 
	calculation of F_prime
Check https://projecteuler.net/problem=625
      https://codeforces.com/problemset/problem/2020/F
*/ const int max_h=1e6;
template<typename tint> struct min25{min25()=default;
vector<int> prm;void get_primes(int n){
    vector<bool> crb(n,0);rep(i,2,n){
        if(!crb[i])prm.pb(i),mu[i]=-1,ph[i]=i-1;
        for(int p:prm){if(p*1ll*i>=n)break;
            crb[i*p]=1;if(i%p==0)break; } }
}int r_up[max_h],r_dwn[max_h];
lint blck[max_h],n;
inline int ren(lint x){return x<max_n?r_dwn[x]:r_up[n/x];}
tint F_prime[max_h],F[max_h]; void solve(){ // make blocks
	int k=0;{lint l=1;while(l<=n){blck[k++]=n/l;l=n/(n/l)+1;}}
    sort(blck,blck+k);rep(i,0,k){
        if(blck[i]<max_n)r_dwn[blck[i]]=i;
        else r_up[n/blck[i]]=i;
    } // for the below line, I need to make $F\_prime[i]=\sum_{j=2}^{blck[i]}g(j)$ for some
    // (potentially other) multiplicative function g. Notice that this misterious g must 
    // be g(p)=f(p), yet in other values it may be completely different. You may actually
    // want to divide this computation in different steps if there is not a single g you can use.
    rep(i,0,k)F_prime[i]=?; // for instance, $\phi(p)=Id(p)-\zeta(p)$
    rep(i,0,prm.size())rep(j,k-1,-1){if(prm[i]*1ll*prm[i]>blck[j])break;
    	F_prime[j]-=g(prm[i])*(F_prime[ren(blck[j]/prm[i])]-(i?F_prime[ren(prm[i-1])]:mint(0)));
	}rep(i,0,k)F[i]=F_prime[i];
	rep(i,prm.size()-1,-1){ const int p=prm[i],rp=ren(p);
        rep(j,k-1,-1){ if(p*1ll*p>blck[j])break;
            for(lint pc=p;pc*p<=blck[j];pc*=p)
                F[j]+=f(pc)*(phi[ren(blck[j]/pc)]-F_prime[rp])+f(pc*p);
        }
	} rep(i,0,k)F[i]+=1;
} };