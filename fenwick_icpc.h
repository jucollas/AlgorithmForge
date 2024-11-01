/*
Autor: Oscar Vargas Pabon
Material de referencia para ICPC
Lo probe en testing\rangeQtest.cpp
*/

int bit[template_limit], bit_size; // BIT -> Binary Indexed Tree

void build( int *arr, int n ) {
	/* Crea el arbol en tiempo O(n) sobre el arreglo 'arr' */
	int nxt; bit_size = n+1; // porque esta indexado en 1
	memset( bit, 0, sizeof(int)*bit_size );
	for ( int i = 0 ; i < bit_size ; ++i ) {
		bit[i] += arr[i]; nxt = i + (i&(-i));
		if ( nxt < bit_size ) bit[nxt] += bit[i];
	}
}
int query( int x ) {
	/* Responde la suma del prefijo de bit[1]+...+bit[x] en tiempo O(lg n)*/
	int res = 0;
	for ( ; x > 0 ; x -= x&(-x) ) res += bit[x];
	return res;
}
int query( int l, int r ) {
	/* responde la suma del rango de bit[l]+...+bit[r] en tiempo O(lg n) */
	return query( r ) - query( l-1 );
}
void update( int x, int val ) {
	/* Modifica el valor de arr[x] en la representacion arborea en O(lg n) */
	for ( ; x < bit_size ; x += x&(-x) ) bit[x] += val;
}

// computa el logaritmo entero de num(la cantidad de bits activos) - esta fun es aparte
int ilog2( int num ) { return sizeof(int)*8 - __builtin_clz( num ) -1; } 
int binlift( int val ) {
	/* Retorna el menor indice 'x' tal que query( x ) >= val en tiempo O(lg n) */
	int x = 0, nxt, sm = 0;
	for ( int exp = ilog2(bit_size) ; exp >= 0 ; --exp ) {
		nxt = x|(1<<exp);
		if ( nxt < bit_size && bit[nxt] < val ) {
			val -= bit[nxt];
			x = nxt;
		}
	}
	return x+1;
}

