/*
Autor: Oscar Vargas Pabon

Esta untested. Notar que no es trivial deshacerse del grafo transpuesto:
	G:=<{1,2,3},{(1,2),(1,3),(2,1)}>
	Si mi toposort organiza [1,2,3] entonces asumira que {1,2,3} es un SCC
		lo cual no es cierto. (idea de IntroToAlgorithm_CormenEtAl Chapter20,
			excercise 20.5-3, fourth edition)
*/
vector<int> topo(const vector<vector<int>> &g){
	int n=g.size();
	vector<int> res;res.reserve(n);
	vector<bool> vis(n,0);
	function<void(int)> taux=[&](int nd){
		vis[nd]=1;
		for(int e:g[nd])if(!vis[e])taux(e);
		res.push_back(nd);
	}; rep(i,0,n)if(!vis[i])taux(i);
	reverse(res.begin(),res.end()); return res;
}
void kosaraju(const vector<vector<int>> &g){
	const int n=g.size();
	
	vector<vector<int>> ig(n);// Transpose graph
	for(int i=0;i<n;++i)for(int e:g[i])ig[e].push_back(i);
	
	vector<int> tord=topo(g);
	vector<int> gr(n,-1);// finding group
	function<void(int)> dfs=[&](int nd,int gind)->void{
		gr[nd]=gind;for(int e:g[nd])if(gr[e]==-1)dfs(e,gind);
	}; int gind=0;for(int nd:tord)if(gr[nd]==-1)dfs(nd,gind++);
	vector<vector<int>> cg(gind);//condensed graph
	for(int i=0;i<n;++i)for(int e:g[nd])
		cg[gr[i]].push_back(gr[e]);
}