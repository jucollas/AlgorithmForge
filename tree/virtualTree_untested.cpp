/*
Autor: Oscar Vargas Pabon
Material de referencia para ICPC
NO LO HE PROBADO

La mayor idea de esto es lo de que para todo $A\subseteq V$ se cumple que $H=\{lca(u,v)|u,v\in A\}$
	es igual a yo ordenar por tiempo e entrada del dfs y solo hacer $\{lca(A_i,A_{i+1})\}_{i=0}^{|A|-1}$

Asumo un lca con funcion 'dfs_time', 'build_lca' y arreglos 'tin', 'tout'(este
		ultimo habria que agregarlo a mano a la implementacion de LCA propuesta)
*/

bool isParent( int u, int v ) {
	/* Responde si 'u' es padre de 'v' en el arbol */
	bool res = tin[u] < tin[v] && tout[u] >= tout[v];
	return res;
}

vector<pair<int,list<int>>> build_virtualTree( const vector<int> &vertex ) {
/* Retorna un arbol virtual en O(s lg s) con raiz en 0 */
	if ( vertex.empty() ) return vector<pair<int,list<int>>>();
	int n = vertex.size();
		
	vector<pair<int,int>> vtNode(n); // para hacer un sort por los tin de los vertices
	for ( int i = 0 ; i < n ; ++i ) vtNode[i] = pair<int,int>( tin[vertex[i]], vertex[i] );
	sort( vtNode.begin(), vtNode.end() );
	
	for ( int i = 1 ; i < n ; ++i ) { // añadir los LCA necesarios para especificar el arbol
		int new_vertex = lca( vtNode[i-1].second, vtNode[i].second );
		vtNode.push_back( pair<int,int>( tin[new_vertex], new_vertex ) );
	}
	sort( vtNode.begin(), vtNode.end() );
	
	int implementationVertex = n; // only to make the implementation easier
	tin[n] = -1; tout[n] = n+1;
	vtNode.push_back( pair<int,int>( -1, implementationVertex ) );
		
	vector<pair<int,list<int>>> tree ( 1, pair<int,list<int>>( vtNode.front().second, list<int>() ) );
	
	map<int,int> rename; // el nombre del vertice en el virtual-tree
	rename[vtNode[0].second] = 0;
	
	list<int> vertex_stack; // el stack para construir el vt
	vertex_stack.push_back( vtNode.front().second );
	
	for ( int i = 1 ; i <= n ; ++i ) { // añade las condiciones 
		if ( vtNode[i].second == vertex_stack.back() ) continue; // to prevent overcounting
		
		// añadir el nuevo vertice
		rename.insert( pair<int,int> ( vtNode[i].second, rename.size() ) );
		tree.push_back( pair<int,list<int>> ( vtNode[i].second, list<int>() ) );
			
		while ( int(vertex_stack.size()) >= 2 && !isparent(vertex_stack.back(),vtNode[i].second) ) {
			int tmp = rename[vertex_stack.back()]; vertex_stack.pop_back();
			
			// añado la arista entre este y el anterior elemento
			tree[tmp].second.push_back( rename[vertex_stack.back()] );
			tree[rename[vertex_stack.back()]].second.push_back( tmp );
		}
		vertex_stack.push_back( vtNode[i].second );
	}
	tree.pop_back(); // para no contar el 'implementationVertex'
	return tree;
}

