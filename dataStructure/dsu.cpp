/*
Autor: Oscar Vargas Pabon
Material de referencia para ICPC
Lo probe en https://codeforces.com/problemset/problem/1857/G
*/

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
	