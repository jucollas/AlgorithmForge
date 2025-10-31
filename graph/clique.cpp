/*
Autor: Oscar Vargas Pabon

Backtracking version tested in https://open.kattis.com/problems/maxclique
	(it seems to do well with V<=50)
Diamond free version tested in https://codeforces.com/gym/105505/problem/D

Add version for chordal graphs
*/

// Note it can be transformed to give the maximal clique by
// keeping the used vertices instead of 'am'
int max_clique(const vector<lint> &g){
	// I'm assuming no self-edges
	int res=0,n=g.size();
	auto ilog2ll=[]( lint num ) { return 8*sizeof(lint) - __builtin_clzll( num ) - 1; };
	function<void(int,lint)> backt=[&](int am,lint msk){
		if(am+__builtin_popcountll(msk)<=res)return;
		
		if( msk==0ll )res=max(res,am);
		else {
			int i=ilog2ll(msk&(-msk));
			backt(am+1,msk&g[i]);
			backt(am,msk&(~(1ll<<i)));
		}
	};backt(0,(1ll<<n)-1);
	return res;
}

// For a diamond free graph, each edge belongs to,
// at most, 1 maximal clique. O(n^2)
vector<vector<int>> chordal_clique(const vector<vector<bool>> &g){
	// I assume adj-matrix representation of g
	int n=g.size();vector<vector<int>> res;
	vector<vector<bool>> vis(n,vector<bool>(n,0));
	rep(i,0,n)rep(j,0,n)if(!vis[i][j] && g[i][j] ){
		vector<int> act={i,j}; //now I search the group
		rep(k,0,n)if(g[i][k]&&g[j][k])act.pb(k);
		for(int u:act)for(int v:act)vis[u][v]=1;
		res.pb(act);
	}
	return res;
}