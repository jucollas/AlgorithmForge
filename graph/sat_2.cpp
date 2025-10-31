/*
Autor: Oscar Vargas Pabon

Usado en https://codeforces.com/contest/2120/problem/F
y en https://codeforces.com/gym/105053/problem/E

Esta version no esta testeada por completo (es mas para dar la idea)
Mi representacion de las clausulas es una mierda
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
bool check_sat(const vector<pair<pair<int,bool>,pair<int,bool>> &claus, const vector<bool> &ass){
	bool res=0;
	for( pair<pair<int,bool>,pair<int,bool>> ac:claus){
		res=res&&((ass[ac.first.first]^(!ac.first.second))||(ass[ac.second.first]^(!ac.second.second)));
	}
	return res;
}
vector<bool> sat2(const vector<pair<pair<int,bool>,pair<int,bool>>> &claus,int n){
	vector<vector<int>> g(n);
	for(pair<pair<int,bool>,pair<int,bool>> ac:claus)rep(i,0,2){
		g[ac.first.first*2+ !ac.first.second].push_back(ac.second.first*2 + ac.second.second);
	}
	vector<int> tord=topo(g);
	vector<int> vl(n,-1);
	rep(i,n*2-1,-1) if ( vl[tord[i]/2]==-1 ) vl[tord[i]/2]=tord[i]&1;
	vector<bool> res(n);rep(i,0,n)res[i]=vl[i];
	return (check_sat(claus,res))?res:{};
}