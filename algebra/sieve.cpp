/* Author: Oscar Vargas Pabon

Based loosely on https://codeforces.com/blog/entry/54090
Can be extended to calculate other multiplicative functions such as
	$\phi(p^e)=p^{e-1}(p-1)$ or ad-hoc stuff (remember I can mantain the value of e
	to calculate certain weird multiplicative functions)
*/
int mu[max_n];void sieve(int n){
	vector<bool>crb(n,0);vector<int>prm;
    crb[1]=mu[1]=1;rep(i,2,n){
        if(!crb[i])prm.pb(i),mu[i]=-1;
        for(int p:prm){if(i*1ll*p>=n)break;
            crb[i*p]=1;
            if(i%p==0){mu[i*p]=0;break;}
            mu[i*p]=mu[i]*mu[p];
        }
    }
}