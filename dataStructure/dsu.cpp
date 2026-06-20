/*
Autor: Oscar Vargas Pabon
Recordar que existe el 'reachibility tree'
	(dejar las uniones explicitas como nuevos nodos, sacrificando compresion de caminos)
TODO: En la superregional se hizo un 'augmenting' para que acumulara en sets.
Lo probe en https://codeforces.com/problemset/problem/1857/G
*/
struct dsu{ // currently untested
	int pi[max_n],sz[max_n],n;
	inline int add_set(){pi[n]=n;sz[n]=1;return n++;}
	dsu(int nn=0){n=0;while(nn--)add_set();}
	int get_pi(int x){return pi[x]=pi[x]==x?x:get_pi(pi[x]);}
	bool merge(int x,int y){
		x=get_pi(x);y=get_pi(y); if(x==y)return 0;
		if(sz[x]>sz[y])swap(x,y);
		pi[x]=y;sz[y]+=sz[x]; return 1; }
};

//////////// previous version

int dsu_pi[template_limit], dsu_sz[template_limit];

void build( int n ) {
/* Construye el estado inicial en tiempo O(n) */
	for ( int i = 0 ; i < n ; ++i ) dsu_pi[i]=i, dsu_sz[i]=1;
}

int getRepr( int x ) {
/* Halla el representante de x */
	if ( dsu_pi[x] != x ) dsu_pi[x] = getRepr( dsu_pi[x] );
	return dsu_pi[x];
}

void mergeSet( int x, int y ) {
/* Une los sets de 'x' y 'y' */
	x = getRepr( x ); y = getRepr( y );
	if ( x != y ) {
		if ( dsu_sz[x] < dsu_sz[y] ) swap( x, y );
		dsu_pi[y] = x;
		dsu_sz[x] += dsu_sz[y];
	}
}
