/*
Autor: Oscar Vargas Pabon
Material de referencia para ICPC
Lo probe en mi maquina
Esta implementacion asume una RMQ-DS con funciones
	* void build( int *arr, int n ) * int query(int l, int r )
*/

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
