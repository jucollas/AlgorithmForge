/*
Autor: Oscar Vargas Pabon
Material de referencia para ICPC
Lo probe con testing\rangeQtest.cpp

lt -> left tree bound, rt -> right tree bound
lx -> left op bound, rx -> right op bound
*/

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
		tag[ind+2*(((rt+lt)/2)-lt+1)] += tag[ind]; // posible issue at optimization
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