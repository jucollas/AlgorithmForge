int pi[max_n][2],dpt[max_n];
void build_lca(const vector<vector<int>>&t,int nd=0,int p=0,int d=0){
	if(!nd)pi[0][0]=pi[0][1]=0;
	dpt[nd]=d;pi[nd][0]=p; // this ensures somehow O(lgn) query time
	int pp=pi[p][1],ppp=pi[pp][1]; if(dpt[p]-dpt[pp]==dpt[pp]-dpt[ppp])
		pi[nd][1]=ppp;
	else pi[nd][1]=p,bl[nd][1]=arr[nd];
	for(int e:t[nd])if(e!=p)build_lca(t,e,nd,d+1);
} int lca(int u,int v){ if(dpt[u]<dpt[v])swap(u,v);
	while(dpt[u]>dpt[v])u=pi[u][dpt[pi[u][1]]>=dpt[v]];
	while(u!=v){ int jmp=pi[u][1]!=pi[v][1];
		u=pi[u][jmp];v=pi[v][jmp];
	} return u;//notice I am already on u
} inline int dist(int u,int v){// counts distance in edges
    const int p=lca(u,v); //+1 to count in nodes
    return dpt[u]+dpt[v]-2*dpt[p]; }