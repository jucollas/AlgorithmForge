/*
Autor: Oscar Vargas Pabon
*/

#include <bits/stdc++.h>
using namespace std;

const int template_limit = 1e6;
/* insertar RMQ-DS */

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

/* Insertar LCA */

int euler[template_limit], tin[template_limit], itin[template_limit];

int dfs_time( const vector<list<int>> &tree, int node, int father=-1, int time=0 ) {
	tin[node] = time; itin[time] = node; // los renombramientos
	euler[time] = time; // el orden en la RMQ-DS
	++time;
	
	for ( const int it : tree[node] ) {
		if ( it != father ) {
			time = dfs_time( tree, it, node, time );
			euler[time] = tin[node]; // para que aparezca entre todos sus sub-arboles
			++time;
		}
	}
	return time;
}
void build_lca( const vector<list<int>> &tree, int root=0 ) {
/* construyo la LCA en tiempo O(n+B(n)) donde B es el tiempo de construir de la RMQ-DS */
	int tmp = dfs_time( tree, root ); // calculo los tiempos y -euler circuit-
	build( euler, tmp ); // construyo la RMQ-DS
}
int lca( int u, int v ) {
/* respondo el LCA en tiempo O(T(n)) donde T es el tiempo de query de la RMQ-DS */
	u = tin[u]; v = tin[v];
	if ( u > v ) swap( u, v );
	return itin[query(u,v)];
}

/* el resto */

int random( int tope ) {
        int i = abs( rand() % tope );
        return i;
}
vector<list<int>> treeMaker( int sz ) {
	vector<list<int>> tree(sz);
	for ( int i = 1 ; i < sz ; ++i ) {
		int pr = random( i );
		tree[pr].push_back( i );
		tree[i].push_back( pr );
	}
	return tree;
}

pair<int,int> lca_brute( const vector<list<int>> &tree, int node, int parent, int a, int b ) {
	pair<int,int> res(-1,0);
	for ( const int edge : tree[node] ) {
		if ( edge != parent ) {
			pair<int,int> tmp = lca_brute( tree, edge,node, a, b );
			if ( tmp.first != -1 ) {
				res.first = tmp.first; 
				break;
			}
			res.second |= tmp.second;
		}
	}
	assert ( node != -1 );
	if ( node == a ) res.second |= 2;
	if ( node == b ) res.second |= 1;
	if ( res.second == 3 && res.first == -1 ) res.first = node;
	return res;
}

const int rep = 1e3, tope = 1e5;

int main(){
	while ( 1 ) {
		vector<list<int>> tree = treeMaker( tope );
		int root = random( tope );
		build_lca( tree, root );
		for ( int i = 0 ; i < rep ; ++i ) {
			int a = random( tope), b = random(tope);
			if ( lca( a, b ) != lca_brute( tree, root, -1, a, b ).first ) {
				cout << "root -> " << root << endl;
				for ( int i = 0 ; i < tope ; ++i  ) {
					cout << "[" << i << "]";
					for ( const int edge : tree[i] ) cout << ' ' << edge ;
					cout << endl;
				}
				cout << "("<<a<<","<<b<<") -> " << lca(a,b) << ' ' << lca_brute(tree,root,-1,a,b).first << endl;
				
			}
			assert( lca(a,b) == lca_brute( tree, root, -1, a, b ).first );
		}
	}
	return 0;
}