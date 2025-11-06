/*
Autor: Oscar Vargas Pabon
Material de referencia para ICPC
Lo probe en https://codeforces.com/problemset/problem/2026/E 
Notar las variables globales:
    * vector<list<int>> graph;
	* vector<Edge> edges;
	* int source, sink;
	*
	* int level[template_limit];
	* bool blocked[template_limit];
	* list<int>::iterator ptr[template_limit];
Mi implementacion interpreta los valores de 'graph' como indices de 'edges'
*/

class Edge {
public:
	int u, v, c1, c2, f;// c1 es la capacidad u->v y c2 es para v->u
	Edge()=default;
	Edge( int u, int v, int c1, int c2=0 ) : u(u),v(v),c1(c1),c2(c2),f(0) {};
	int to( int node ) const { // el opuesto a 'node'
		return ( u == node ) ? v : u;
	}
	void push( int node, int pushed ) { // empujar un flujo 'pushed' desde 'node'
		if ( node == u ) c1 -= pushed, c2 += pushed, f += pushed;
		else c1 += pushed, c2 -= pushed, f -= pushed;
	}
	int cap( int node ) const { // la capacidad desde 'node'
		return ( this->u == node ) ? c1 : c2;
	}
	int flow( int node ) const { // el flujo desde 'node'
		return ( u == node ) ? f : -f;
	}
};

int level[template_limit]; // para el 'layered network'
bool blocked[template_limit]; // para omitir vertices 'blocked'
list<int>::iterator ptr[template_limit]; // para omitir aristas 'blocked'

vector<list<int>> graph;
vector<Edge> edges;
int source, sink;

void addEdge( int u, int v, int c1, int c2=0 ) {
	// u->v has capacity c1; v->u has capacity c2
	
	graph[u].push_back( edges.size() );
	graph[v].push_back( edges.size() );
	edges.push_back( Edge( u, v, c1, c2 ) );
}

bool level_bfs( ) {
/* construye el 'layered network' de manera implicita con 'level' en O(V+E) */
	memset( level, -1, sizeof(int)*int(graph.size()) );
	level[source] = 0; // 'level' es la profundidad el el BFS-Tree
	
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
	}
	return level[sink]!=-1; // para saber si ya terminamos
}

int push_dfs( int node, int flow=1e9 ) {
/* Empuja el flujo por todas las aristas del 'layered graph' que pueda */
	if ( node == sink || flow==0 ) return flow; // ya llegamos o no podemos empujar mas
	
	int push_flow = 0, edge_flow, edge, nxt;
	while ( ptr[node] != graph[node].end() && flow > push_flow ) {
		edge_flow = 0; edge = *ptr[node]; nxt = edges[edge].to(node);
		// si la arista pertenece al 'layered graph' y el vertice NO esta bloqueado
		if ( level[node] < level[nxt] && !blocked[nxt] )  {
			edge_flow = push_dfs( nxt, min(edges[edge].cap(node),flow-push_flow) );
			push_flow += edge_flow;
			edges[edge].push( node, edge_flow );
		}
		++ptr[node];
	}
	--ptr[node]; // para corregir el ultimo ++ptr[node]
	blocked[node] = (push_flow == 0); // verifica si se bloqueo el vertice
	return push_flow;
}

int maxFlow( ) {
/* Hace el maximo flujo del grafo (retorna el maxFlow, las asignaciones quedan en las aristas)
	Funciona en peor caso $O(v^2*E)$ */
	int flow = 0;
	while ( level_bfs( ) ) {
		memset( blocked, 0, sizeof(bool)*int(graph.size()) ); // reinicializo esto
		for ( int i = 0 ; i < int(graph.size()) ; ++i ) ptr[i] = graph[i].begin();
		
		flow += push_dfs( source );
	}
	return flow;
}