/*
Autor : Oscar Vargas Pabon
Material de referencia para ICPC
Lo probe en 12299 - RMQ with Shifts
*/
int tree[template_limit], arr_size;

void build( int *arr, int n ) {
	/* Construye el arbol en tiempo O(n) */
	arr_size = n;
	for ( int i = 0 ; i < arr_size ; ++i ) tree[i+n] = arr[i];
	for ( int i = arr_size-1 ; i > 0 ; --i ) tree[i] = min( tree[i*2], tree[i*2+1] );
}
void update( int x, int val ) {
	/* Modifica el valor en arr[x] segun la representacion del arbol en O(lg n) */
	x+=arr_size;
	tree[x]=val;
	for ( x >>= 1 ; x > 0 ; x >>= 1 ) tree[x] = min( tree[x*2], tree[x*2+1] );
}

int query( int l, int r ) {
	/* Responde a la query arr[l] + ... + arr[r] en tiempo O(lg n) */
	int res = 1e9 ;
	for ( l += arr_size, r += arr_size ; l <= r ; l >>= 1, r >>= 1 ) {
		if ( l&1 ) res = min( res, tree[l++] );
		if ( !(r&1) ) res = min( res, tree[r--] );
	}
	return res;
}