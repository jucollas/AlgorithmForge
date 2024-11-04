/*
Autor: Oscar Vargas Pabon
Material de referencia para ICPC

Probado en CSES 2134 -- 1138

Esta implementacion asume una Range-DS con funciones
	* void build( int *arr, int n ) * int query(int l, int r )
*/

int tin[template_limit], tout[template_limit], sz[template_limit]; // data del arbol
int jmp[template_limit], chain[template_limit]; // data de los saltos de HLD
int perm[template_limit]; // data para la Range-DS

bool isParent( int u, int v ) {
	/* Responde si 'u' es padre de 'v' en el arbol */
	bool res = tin[u] < tin[v] && tout[u] >= tout[v];
	return res;
}

void getSz( vector<list<int>> &tree, int node=0, int parent=-1 ) {
	/* Hallo el tamaño, pongo las Heavy-edges de primero */
	sz[node] = 1;
	int heavy = node;
	for ( const int it : tree[node] ) {
		if ( it != parent ) {
			getSz( tree, it, node );
			if ( sz[heavy] < sz[it] ) heavy = it; // Heavier edge
			sz[node] += sz[it]; // calculo el tamaño
		}
	}
	if ( heavy != node ) { // reorganizo para que las Heavy-edges queden de primero
		tree[node].remove( heavy );
		tree[node].push_front( heavy );
	}
}

void getLabeling( const vector<list<int>> &tree, int node=0, int parent=-1, int time=0 ) {
	/* Hallo los tin,tout, y los jmp,chain */
	tin[node] = time++; // hallo el tiempo de entrada
	bool first = true;
	for ( const int it : tree[node] ) {
		if ( it != parent ) {
			if ( first ) { // esta es una Heavy-Edge
				jmp[it] = jmp[node];
				chain[it] = chain[node];
				first = false;
			} else { // las demas
				jmp[it] = node;
				chain[it] = it;
			}
			
			getLabeling( tree, it, node, time );
			time = tout[it]; // actualizo el tiempo
		}
	}
	tout[node] = time; // encuentro el tout
}

void build_hld( vector<list<int>> &tree, const vector<int> &value ){
	/* Construye el hld a partir de un arbol y sus valores asociados en tiempo
		O(n+B(n)) donde B indica el tiempo de la Range-DS usada en inicializar */
	const int root = 0; // arbitrario
	getSz( tree, root ); // hallo las Heavy-Edges
	jmp[root] = root; 
	getLabeling( tree, root ); // hallo tin,tout,jmp y chain
	
	// construyo la Range-DS
    for ( int i = 0 ; i < n ; ++i ) perm[tin[i]] = value[i];
    build( perm, n );
}

int queryPath( int u, int v ) {
	/* Responde a la query arr[u] + .. + arr[v] donde todos los elementos estan en el 
		camino de 'u' a 'v' en el arbol. En tiempo O(lg n *T(n) ) donde T(n) es el 
		tiempo de la DS usada */
	int res = 0;
	while ( chain[u] != chain[v] ) {
		if ( isParent( chain[u], v ) ) swap( u, v ); // el que llegue primero hace swap
		res += query( tin[chain[u]], tin[u] );
		u = jmp[u];
	}
	
	if ( isParent( u, v ) ) swap( u, v ); // para incluir el chain hasta el LCA
	res += query( tin[v], tin[u] );
	return res;
}