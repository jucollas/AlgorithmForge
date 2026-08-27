/* Autor: Oscar Vargas Pabon
Lo probe en https://codeforces.com/problemset/problem/2026/E 
Mi implementacion interpreta los valores de 'graph' como indices de 'edges'
Caso general funciona en $O(v^2*E)$, caso capacidad unitaria funciona en 
	$O(min(E\sqrt{E},EV^{2/3}))$
	caso red unitaria (cada nodo se responsabiliza por 1 sola arista de
		capacidad unitaria) funciona en $O(E\sqrt{V})$

Remember,remember. Whenever we have a planar graph, a path of its dual graph represents
	bijectively a cut on the original graph. Essentially solves with dijkstra min-cut
Anfelesan taught me this one in regards to https://codeforces.com/gym/106178/problem/G
*/ template <typename ftype, ftype INF> struct Dinitz{
	constexpr static ftype inf(){return INF;}
	struct Edge { int u, v;
		ftype c1, c2, f;// c1 es la capacidad u->v y c2 es para v->u
		Edge()=default;
		Edge( int uu, int vv, ftype cc1, ftype cc2=0 ) : u(uu),v(vv),c1(cc1),c2(cc2),f(0) {};
		inline void push( int nd, ftype pushed ) { // empujar un flujo 'pushed' desde 'nd'
			if ( nd == u ) c1 -= pushed, c2 += pushed, f += pushed;
			else c1 += pushed, c2 -= pushed, f -= pushed;
		} inline int to( int nd ) const { // el opuesto a 'nd'
			return ( u == nd ) ? v : u;
		} inline ftype cap( int nd ) const { // la capacidad desde 'nd'
			return ( this->u == nd ) ? c1 : c2;
		} inline ftype flow( int nd ) const { // el flujo desde 'nd'
			return ( u == nd ) ? f : -f; }
	}; // level -> para el 'layered network'
	// blocked-> para omitir vertices 'blocked'
	// ptr -> para omitir aristas 'blocked'
	vector<int> level,ptr; vector<bool>blocked;

	vector<vector<int>> graph; vector<Edge> edges;
	int source, sink; Dinitz(int n){
		graph.resize(n); level.resize(n);
		ptr.resize(n); blocked.resize(n);
	} void add_edge( int u, int v, ftype c1, ftype c2=0 ) {
		// u->v has capacity c1; v->u has capacity c2	
		graph[u].push_back( edges.size() );
		graph[v].push_back( edges.size() );
		edges.emplace_back( Edge( u, v, c1, c2 ) );
	} bool level_bfs( ) {
// construye el 'layered network' de manera implicita con 'level' en O(V+E) 
		fill(all(level),-1); level[source] = 0;
		// 'level' es la profundidad el el BFS-Tree
		queue<int> q; q.push( source );
		while ( !q.empty() && level[sink] == -1 ) {
			int act = q.front(); q.pop();
			for ( int edge : graph[act] ) {
				int nxt = edges[edge].to(act); // el siguiente vertice
				if ( level[nxt] == -1 && edges[edge].cap(act)>0 ) {
					level[nxt] = level[act] + 1;
					q.push( nxt );
				}
			}
		} return level[sink]!=-1; // para saber si ya terminamos
	} ftype push_dfs( int nd, ftype flow=INF ) {
// Empuja el flujo por todas las aristas del 'layered graph' que pueda 
		if ( nd == sink || flow==0 ) return flow; // ya llegamos o no podemos empujar mas
		
		ftype push_flow = 0, edge_flow; int edge, nxt;
		while ( ptr[nd] < int(graph[nd].size()) && flow > push_flow ) {
			edge_flow = 0; edge = graph[nd][ptr[nd]]; nxt = edges[edge].to(nd);
			// si la arista pertenece al 'layered graph' y el vertice NO esta bloqueado
			if ( level[nd] < level[nxt] && !blocked[nxt] )  {
				edge_flow = push_dfs( nxt, min(edges[edge].cap(nd),flow-push_flow) );
				push_flow += edge_flow;
				edges[edge].push( nd, edge_flow );
			} ++ptr[nd];
		} --ptr[nd]; // para corregir el ultimo ++ptr[nd]
		blocked[nd] = (push_flow == 0); // verifica si se bloqueo el vertice
		return push_flow;
	} ftype max_flow( ) {
// Hace el maximo flujo del grafo (retorna el maxFlow, las asignaciones quedan en las aristas)
//		Funciona en peor caso $O(V^2*E)$
		ftype flow = 0; while ( level_bfs( ) ) {
			fill(all(blocked),0); fill(all(ptr),0);
			
			flow += push_dfs( source );
		} return flow; }
}; typedef Dinitz<lint,lint(1e18+333)> Flow; typedef Flow::Edge Edge;