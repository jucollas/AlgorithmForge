/* Autor:Oscar Vargas Pabon
DSU que permite hacer rollback de la ultima operacion (como un stack)
Funciona en tiempo O(lg n) porque no se puede hacer 'path-compression'
Util en cosas como offline dynamic connectivity
Probado con impl de offlineDynamicConnectivity en
	https://codeforces.com/gym/100551/problem/A
Me gusto un problema de offlineDynamicConnectivity en
	https://www.luogu.com.cn/problem/P5625
Remember that to do XOR path queries (given weighted edges) I do
	int path_sum(int x){int x=0;if(pi[x]!=x)x=wg[x]^path_sum(pi[x]);return x;}
	And in merge I modify 
		wg[x]=edge_weight^path_sum(ex)^path_sum(ey); pi[x]=y; cnt[y]+=cnt[x];
	where (ex,ey) is the added edge
Can I extend it to + path queries or * path queries????
*/ struct DSU{ vector<int> pi,cnt,stc;
	DSU(int n=0){ pi.resize(n);cnt.resize(n);stc.resize(0);
		for(int i=0;i<n;++i)pi[i]=i,cnt[i]=1;
	} int get_pi(int nd){return pi[nd]==nd?nd:get_pi(pi[nd]);}
	bool merge(int x,int y){
		x=get_pi(x);y=get_pi(y);
		if(x!=y){ if(cnt[x]>cnt[y])swap(x,y);
			pi[x]=y; cnt[y]+=cnt[x];
			stc.push_back(x);
		}else stc.push_back(-1);
		return x!=y;
	} void rollback(int am=1){
		while(am--){ assert(!stc.empty());
			int x=stc.back(); stc.pob();
			if(x!=-1){ cnt[pi[x]]-=cnt[x];
				pi[x]=x; }
		}
	} };