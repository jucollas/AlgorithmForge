/*
Author: Oscar Vargas Pabon
ChineseRemainderTheorem

Helps to represent numbers muli composites with $O(lg m)$ operations.
I added a weird impl that generates nChoose k in time 
	$O(log_p(n))$ for each prime p. This uses the idea of generalized lucas theorem
	
This was based on https://www.luogu.com.cn/article/g3wlhtkg
	* https://www.codeleading.com/article/4340352347/
	* https://www.cnblogs.com/JohnYam/p/19055958
Tested on https://www.luogu.com.cn/problem/P3301

Note that I could simplify a bit the CRT having something like 
ll crt(ll a,ll pk) { ll x=p/pk; return a*x%p*inv(x,pk)%p; }
invoked once per prime component of the modulus

My impl depends on 
crt::comb  ->  mpow<tpow>, crt_info.fact
All -> crt_info.pk , crt_info.p

Note that with a little more of effort, this could become an arbitrary modulus template.
*/
namespace internal{

template <typename tgcd>
tgcd exgcd (tgcd a, tgcd b, tgcd &x, tgcd &y) {
	if (b==0) { x = 1;y = 0; return a; }
	tgcd d = exgcd(b, a % b, x, y), t = x;
	x = y, y = t - a / b * y;
	return d;
}
template <typename tinv>
inline tinv exinv(tinv x,tinv m){
	// $x*s+m*t=g$
	tinv s,t,g;g=exgcd<tinv>(x,m,s,t); (void)g; (void)t;
	return ((s%m)+m)%m;
}
inline int padic_fact(lint x,lint p){
	//$res = v_p(x!)$ 
	//using legendre's formula $res=\sum_{e>0}\lfloor n/(p^k)\rfloor$
	int res=0;lint pp=1;while(x>pp){
		pp*=p;res+=x/pp;
	}
	return res;
}

template<typename tcrt, tcrt m>
struct crt_info{ static_assert(m>0);
	int n;
	vector<tcrt> p,pk,k;// p->prime; p^k->mulo ; k->degree
	vector<vector<tcrt>> fact;
	crt_info(){
		tcrt md=m; // factorization part
		for(tcrt ind=2;ind*ind<=md;++ind)if(md%ind==0){
			p.push_back(ind);k.push_back(0);pk.push_back(1);
			while(md%ind==0)pk.back()*=ind,++k.back(),md/=ind;
		} if(md>1)p.push_back(md),k.push_back(1),pk.push_back(md);
		n=p.size();
		
		for(int ac=0;ac<n;++ac){ // precalculating factorials
			fact.push_back(vector<tcrt>(pk[ac]));
			fact[ac][0]=1;for(tcrt i=1;i<pk[ac];++i)fact[ac][i]=(i%p[ac])?(fact[ac][i-1]*1ll*i)%pk[ac]:fact[ac][i-1];
		}
	}
};

template<typename tcrt,tcrt m>
struct crt{
	static constexpr tcrt mod(){return m;}
	static const crt_info<tcrt,m> info;
	vector<tcrt> x;
	crt(){x.resize(info.n,0);};
	crt(tcrt v):crt(){v%=m;if(v<0)v+=m; for(int i=0;i<info.n;++i)x[i]=v%info.pk[i]; }
	
	tcrt compute_crt()const{
		tcrt res=0;for(int i=0;i<info.n;++i){
			tcrt mul=m/info.pk[i];
			res+= exinv<tcrt>(mul,info.pk[i])*1ll*mul %m * 1ll * x[i] % m;
			if(res>=m)res-=m;
		}
		return res;
	}
	
	crt &operator *=(const crt&ot){for(int i=0;i<info.n;++i)x[i]=(x[i]*1ll*ot.x[i])%info.pk[i]; return *this;}
	crt  operator * (const crt&ot)const{return crt(*this)*=ot;}
	crt &operator +=(const crt&ot){for(int i=0;i<info.n;++i){x[i]+=ot.x[i];if(x[i]>=info.pk[i])x[i]-=info.pk[i];} return *this;}
	crt  operator + (const crt&ot)const{return crt(*this)+=ot;}
	crt &operator -=(const crt&ot){for(int i=0;i<info.n;++i){x[i]-=ot.x[i];if(x[i]<0)x[i]+=info.pk[i];} return *this;}
	crt  operator - (const crt&ot)const{return crt(*this)-=ot;}
	
	crt operator - ()const{crt res;for(int i=0;i<info.n;++i)res.x[i]=info.pk[i]-x[i];return res;}
	
	bool operator==(const crt&ot)const{return x==ot.x;}
	bool operator!=(const crt&ot)const{return x!=ot.x;}
	bool operator==(const tcrt&ot)const{return x==crt(ot)}
	
	operator int ()const{return compute_crt();}
	operator bool()const{return int(*this);}
	
	
	static crt comb(lint n,lint k){
		// Time $O(\sum_p ( log_p(n) + log(p^k) ) )$
		if(n<k||n<0||k<0)return crt(0);
		
		vector<tcrt> &pk=info.pk,&p=info.p;
		function<tcrt(lint,int)> l_fact=[&](lint v,int ind)->tcrt{
			//extended lucas factorial
			vector<tcrt> &fc=info.fact[ind];
			if(v<p[ind])return fc[v];
			return mpow<tcrt>(fc[pk[ind]-1],v/pk[ind],pk[ind])*1ll*fc[v%pk[ind]]%pk[ind]*1ll
					*l_fact(v/p[ind],ind)%pk[ind];//lucas part -_-
		};
		crt res; for(int i=0;i<info.n;++i){
			int e=padic_fact(n,p[i])-padic_fact(k,p[i])-padic_fact(n-k,p[i]);
			if(e>=info.k[i])continue; assert(e>=0);
			res.x[i]= l_fact(n,i)*1ll*exinv<tcrt>( l_fact(k,i)*1ll*l_fact(n-k,i)%pk[i] ,pk[i] )%pk[i]
						*1ll*mpow<tcrt>(p[i],e,pk[i])%pk[i];
		}
		return res;
	}	
};
}typedef internal::crt<int,mod> crt;