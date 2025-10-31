/*
Autor: Oscar Vargas Pabon
Fecha: 30/10/2024
testeando las impl.
*/
#include <bits/stdc++.h>
using namespace std;

const int template_limit = 1e6;

/* insert here */


int tree[template_limit], tag[template_limit], arr_size;

void raw_build( int *arr, int ind, int lt, int rt ) {
/* funcion recursiva auxiliar que construye el segTree */
	tag[ind] = 0; // siempre en neutro
	if ( lt == rt ) tree[ind] = arr[lt]; // caso base
	else {
		int m = (lt+rt)/2;
		raw_build( arr, ind+1, lt, m ); // hijo izquierdo
		raw_build( arr, ind+2*(m-lt+1), m+1, rt ); // hijo derecho
		
		tree[ind] = tree[ind+1] + tree[ind+2*(m-lt+1)]; // computo este acumulado
	}
}

void build( int *arr, int n ) {
/* construye el segTree de manera recursiva a partir de un arreglo en tiempo O(n) */ 
	arr_size = n;
	raw_build( arr, 0, 0, n-1 );
}

void push( int ind, int lt, int rt ) {
/* Empuja el valor del tag a los nodos mas profundos */
	tree[ind] += tag[ind] * (rt-lt+1); // esto cambia dependiendo de la funcion
	if ( lt < rt ) { // no hay edge-case cuando lt==rt
		tag[ind+1] += tag[ind];
		tag[ind+2*(((rt+lt)/2)-lt+1)] += tag[ind]; // fast-math problem???
	}
	tag[ind] = 0; // para no sobrecontar
}

void update( int lx, int rx, int val, int ind=0, int lt=0, int rt=arr_size-1 ) {
/* actualizo los valores en arr[lx],...,arr[rx] por arr[lx]+=val,...,arr[rx]+=val
	en el arbol en tiempo O(lg n) */
	if ( lx <= lt && rt <= rx ) {
		tag[ind] += val;
		push( ind, lt, rt );
	} else if ( rx >= lt && rt >= lx ) { // se situa en el rango
		push( ind, lt, rt );
		int m = (lt+rt)/2; // el punto medio a cortar
		
		update( lx, rx, val, ind+1, lt, m ); // hijo izquierdo
		update( lx, rx, val, ind+2*(m-lt+1), m+1, rt ); // hijo derecho
		tree[ind] = tree[ind+1] + tree[ind+2*(m-lt+1)]; // actualizo el arbol
	}
}

int query( int lx, int rx, int ind=0, int lt=0, int rt=arr_size-1 ) {
/* Hago una query de arr[lx]+...+arr[rx] en tiempo O(lg n) */
	int res = 0;
	push( ind, lt, rt );
	if ( lx <= lt && rt <= rx ) res += tree[ind];
	else if ( rx >= lt && rt >= lx ) {
		int m = (lt+rt)/2; // el punto medio a cortar
		
		res += query( lx, rx, ind+1, lt, m );
		res += query( lx, rx, ind+2*(m-lt+1), m+1, rt );
	}
	return res;
}

/* test implementation */

const int n = 1e5 * 4, INF = 1e3, indexation = 0;

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

void manualUpdate( int l, int r, int val ) {
	for ( int i = l ; i <= r ; ++ i ) arr[i] += val;
}


int main(){
	
	for ( int i = indexation ; i < n+indexation ; ++i ) arr[i] = random( INF );
	
	build( arr, n );
	for ( int i = indexation ; i < n+indexation ; ++i ) assert ( query(i,i) == arr[i] );
	while ( 1 ) {
		int op = random(2);
		if ( op == 0 ) { // 0 - update
			int x = random(n)+indexation, y=random(n)+indexation, value = random( INF ) ;
			
			value = (arr[x]+value)%INF - arr[x];
			
			update( x, y, value ) ;
			manualUpdate( x, y, value );
			for ( int i = 0 ; i < n ; ++i ) {
				int tmp = query(i,i);
				if( tmp != arr[i] ) {
					cout << "("<<x << ", " << y << " ) -> " << value << endl;
					cout << i << " -> "<< tmp << ", " << arr[i] << endl;
				}
				assert ( query(i,i) == arr[i] );
			}
		} else if ( op == 1 ) { // 1 - query
			int l = random(n)+indexation, r = random(n)+indexation ;
			if ( l> r ) swap( l, r );
			int a = query( l, r ), b = manualQuery( l, r );
			if ( a != b ) {
				cout << "(" << l << ", " << r << ") - " << a << " <-> " << b << endl;
				
			}
			assert( query( l, r ) == manualQuery( l, r ) );
		}
	}
	return 0;
}