/*
Autor: Oscar Vargas Pabon
Material de referencia para ICPC
Lo probe en mi carpeta de pruebas
*/

vector<int> zFunction( const string &cad ) {
/* Computa la funcion z en O(n). 
Esta es z[i] -> maximo prefijo comun de cad y cad[i...] */
	int n = cad.size();
	vector<int> z( n, 0 );
	int l = -1, r = -1;
	for ( int i = 1 ; i < n ; ++i ) {
		z[i] = max( 0, min( r-i, z[i-l] ) );
		while ( i+z[i] < n && cad[z[i]] == cad[i+z[i]] ) ++z[i];
		if ( i+z[i] > r ) r=i+z[i],l=i;
	}
	return z;
}