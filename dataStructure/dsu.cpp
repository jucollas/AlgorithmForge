/* Autor: Oscar Vargas Pabon
Recordar que existe el 'reachibility tree'
	(dejar las uniones explicitas como nuevos nodos, sacrificando compresion de caminos)
TODO: En la superregional se hizo un 'augmenting' para que acumulara en sets.
*/
struct DSU{ vector<int>pi,sz;
	inline int add_set(){pi.push_back(sz.size());sz.push_back(1);return pi.size()-1;}
	DSU(int nn=0,int res=0){res=max<int>(nn,res);pi.reserve(res);sz.reserve(res);while(nn--)add_set();}
	int get_pi(int x){return pi[x]=pi[x]==x?x:get_pi(pi[x]);}
	bool merge(int x,int y){
		x=get_pi(x);y=get_pi(y); if(x==y)return 0;
		if(sz[x]>sz[y])swap(x,y);
		pi[x]=y;sz[y]+=sz[x]; return 1; }
};