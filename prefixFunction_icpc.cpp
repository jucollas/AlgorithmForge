/*
Autor: Oscar Vargas Pabon
Material de referencia para ICPC
Probado en mis biblioteca
*/

vector<int> prefixFunction( const string &cad ) {
/* Computa la prefix function en O(n). Esta es pi[i] -> el tamaño del mayor prefijo 
	de cad que tambien es sujifo de cad[0..i] */
	int n = cad.size();
	vector<int> pi( n, 0 );
	for ( int i = 1 ; i < n ; ++i ) {
		pi[i] = pi[i-1];
		while ( pi[i] > 0 && cad[i] != cad[pi[i]] ) pi[i] = pi[pi[i]-1];
		if ( cad[i] == cad[pi[i]] ) ++pi[i];
	}
	return pi;
}