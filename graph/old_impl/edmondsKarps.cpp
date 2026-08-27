/*
Autor: Oscar Vargas Pabon
Material de referencia para ICPC
Lo probe en https://codeforces.com/problemset/problem/2026/E 
Notar las variables globales:
    * vector<list<int>> graph;
	* vector<Edge> edges;
	* int source, sink;
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
		return ( u == node ) ? c1 : c2;
	}
	int flow( int node ) const { // el flujo desde 'node'
		return ( u == node ) ? f : -f;
	}
};

vector<list<int>> graph;
vector<Edge> edges;
int source, sink;

void addEdge( int u, int v, int c1, int c2=0 ) {
	graph[u].push_back( edges.size() );
	graph[v].push_back( edges.size() );
	edges.push_back( Edge( u, v, c1, c2 ) );
}

int findAugmentingPath(  ) {
/* Hallo un augmenting-path en tiempo O(V+E) y lo actualizo-push */
	vector<int> pred( graph.size(), -1 ); // para hacer una arbol del BFS
	pred[source]=-2; //           Nota: guardamos son las aristas del arbol....
	
	queue<int> q; q.push( source );
	
	while ( !q.empty() && pred[sink] == -1 ) {
		int act = q.front(); q.pop();
		for ( int edge : graph[act] ) {
			int nxt = edges[edge].to( act );
			if ( pred[nxt] == -1 && edges[edge].cap( act ) > 0 ) {
				pred[nxt] = edge; // actualizo la arista usada
				q.push( nxt );
			}
		}
	}
	int flow=-1; // codigo de error
	if ( pred[sink] != -1 ) {
		int act = sink,nxt; flow = 1e9;
		while ( act != source ) { // hallo el flujo
			nxt = edges[pred[act]].to( act );
			flow = min( flow, edges[pred[act]].cap( nxt ) );
			act = nxt;
		}
		act = sink;
		while ( act != source ) { // actualizo las aristas
			nxt = edges[pred[act]].to( act );
			edges[pred[act]].push( nxt, flow );
			act = nxt;
		}
	}
	return flow;
}

int maxFlow(  ) {
/* Calculo el maximo flujo del grafo (retorna el maxFlow, en las aristas queda el flujo) 
	Funciona en peor caso $O(V*(E^2))$ */
	int flow = 0, pFlow=0;
	do {
		flow += pFlow;
		pFlow = findAugmentingPath( );
	} while ( pFlow != -1 );
	return flow;
}