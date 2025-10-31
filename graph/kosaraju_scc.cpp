/*
Autor: Oscar Vargas Pabon

Esta untested. Tengo la gran sospecha de que no es necesario el grafo
	transpuesto, pero no tengo ganas de revisarlo/probarlo
*/
vector<int> topo(const vector<vector<int>> &g){
	int n=g.size();
	vector<int> res; vector<bool> vis(n,0);
	function<void(int)> taux=[&](int nd){
		vis[nd]=1;
		for(int e:g[nd])if(!vis[e])taux(e);
		res.push_back(nd);
	};
	rep(i,0,n)if(!vis[i])taux(i);
	return res;
}
void kosaraju(const vector<vector<int>> &g){
	int n=g.size();
	
	vector<vector<int>> ig(n);// Transpose graph
	for(int i=0;i<n;++i)for(int e:g[i])ig[e].push_back(i);
	
	vector<int> tord=topo(g);
	vector<int> gr(n,-1);int gind=0;// finding group
	function<void(int)> dfs=[&](int nd){
		gr[nd]=gind;for(int e:g[nd])if(gr[e]==-1)dfs(e);
	};
	for(int nd:tord)if(gr[nd]==-1){
		dfs(nd); ++gind;
	}
	vector<vector<int>> cg(gind);//condensed graph
	for(int i=0;i<n;++i)for(int e:g[nd]){
		cg[gr[i]].push_back(gr[e]);
	}
}