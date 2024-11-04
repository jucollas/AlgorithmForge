/*
Autor: Oscar Vargas Pabon
Implementacion de referencia para ICPC

Nomas es una referencia sobre hallar un centroide.
	* Recordar como se puede hacer una 'centroid decomposition' que permite
		hallar propiedades sobre los caminos que pasan por cada centroide asegurando
		O(n lg n )

Probado en CSES 2079
*/

int size[template_limit];

int findSize( const vector<list<int>> &tree, int node=0, int parent=-1 ) {
	/* halla el tamaño de cada subarbol */
	size[node] = 1;
	for ( const int edge : tree[node] ) {
		if ( edge != parent ) size[node] += findSize( tree, edge, node );
	}
	return size[node];
}

int findCentroid( const vector<list<int>> &tree, int node=0 ) {
	/* halla el centroide */
	int centroid = node;
	for ( const int edge : tree[node] ) {
		if ( size[edge] > size[node]/2 ) {
			// reroot to the other edge
			size[node] -= size[edge]; size[edge] += size[node];
			
			centroid = findCentroid( tree, edge );
			break;
		}
	}
	return centroid;
}
