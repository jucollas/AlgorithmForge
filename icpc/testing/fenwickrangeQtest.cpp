/*
Autor: Oscar Vargas Pabon
Fecha: 30/10/2024
testeando las impl.
*/
#include <bits/stdc++.h>
using namespace std;

const int template_limit = 1e6;

/* insert here */

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
	/* Retorna el mayor indice 'x' tal que query( x ) < val en tiempo O(lg n) */
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

/* test implementation */

const int n = 1e5, INF = 1e3, indexation = 1;

int arr[template_limit];

int random( int tope ) {
	int i = abs( rand() % tope );
	return i;
}

int manualQuery( int l, int r ) {
	int res = 0;
	for ( int i = l ; i <= r ; ++i ) res += arr[i];
	return res;
}

void manualUpdate( int x, int val ) {
	arr[x] += val;
}


int main(){
	
	for ( int i = indexation ; i < n+indexation ; ++i ) arr[i] = random( INF );
	int tot = 0;
	for ( int i = indexation ; i < n+indexation ; ++i ) tot += arr[i];
	
	build( arr, n );
	while ( 1 ) {
		int op = random(3);
		if ( op == 0 ) { // 0 - update
			int x = random(n)+indexation, value = random( INF ) ;
			
			value = (arr[x]+value)%INF - arr[x];
			
			update( x, value ) ;
			manualUpdate( x, value );
			tot += value;
			assert( arr[x] >= 0 );
		} else if ( op == 1 ) { // 1 - query
			int l = random(n)+indexation, r = random(n)+indexation ;
			if ( l> r ) swap( l, r );
			// int a = query( l, r ), b = manualQuery( l, r );
			// cout << "(" << l << ", " << r << ") - " << a << " <-> " << b << endl;
			assert( query( l, r ) == manualQuery( l, r ) );
		} else { // 2 - binlift
			int value = random( tot )+1;
			int ind = binlift(value);
			if (!(query( ind-1 ) < value && query(ind)>=value) ) {
				for ( int i = indexation ; i < n+indexation ; ++i ) assert( arr[i]>=0 );
				cout << value << " -> " << ind << endl;
				cout << query( ind - 1 ) << " | " << query(ind) << endl;
			}
			assert( query( ind-1 ) < value && query(ind)>=value );
		}
	}
	return 0;
}