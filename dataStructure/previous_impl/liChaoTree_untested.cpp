/* 
Autor: Oscar Vargas Pabon
Material de referencia para ICPC

NO LO HE PROBADO

CHT
*/

class Line {
	public:
	int m, c;
	Line ( int m, int c ) : m(m),c(c) {};
	int operator() ( int x ) const {
		return m*x+c;
	}
};

Line tree[template_limit];
int tree_size;

void build( int n ) {
	tree_size = n; // inicio todo
	for ( int i = 0 ; i < 2*n ; ++i ) tree[i] = Line( 0, 1e9 );
}

void update( Line &ln, int ind=0, int l=0, int r=tree_size-1 ) {
	int m = ( r+l ) / 2;
	if ( ln(m) < tree[ind](m) ) swap( ln, tree[ind] );
	if ( l < r ) {
		if ( ln(l) < tree[ind](l) ) {
			update( ln, ind+1, l, m );
		} else update( ln, ind+2*(m-l+1), m+1, r );
	}
}

int query( int x, int ind=0, int l=0, int r=tree_size-1 ) {
	/* Recorre el arbol hasta las hojas buscando la mejor raya */
	int m=(r+l)/2;
	int res = tree[ind]( x );
	if ( x < m ) {
		res = min( res, query( x, ind+1, l, m ) );
	} else if ( x > m ) {
		res = min( res, query( x, ind+2*(m-l+1), m+1, r ) );
	}
	return res;
}
