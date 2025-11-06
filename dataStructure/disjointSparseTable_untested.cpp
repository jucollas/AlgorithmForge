/*
Autor: Oscar Vargas Pabon
Fecha: 4/02/2024

|----------------------------------------|
|			  complejidad                |
|----------------------------------------|
|	memoria       |     O( n lg n )      |
|-----------------+----------------------|
|	constructor   |     O( n lg n )      |
|-----------------+----------------------|
|	query         |        O( 1 )        |
|----------------------------------------|

Soporta funciones con propiedad distributiva.
Distributiva: a + ( b + c ) = ( a + b ) + c
*/
// int ilog2( int num ) { return (sizeof(int) * 8) - __builtin_clz( num ) - 1;}
// asumo el anterior de mi template


typedef int Typ;
class DisjointSparseTable {
private:
	static int min( int a, int b ) { return ( a > b ) ? b : a; }
public:
	vector<vector<Typ>> table;
	DisjointSparseTable( const vector<Typ> &vec ) {
		int lg = ilog2( vec.size() ), i, strt, md, j;
		table = vector<vector<Typ>> ( lg+1, vector<Typ> ( vec.size() )  );
		Typ acum;
		for ( i = 0; i < int(vec.size()) ; ++i ) table[0][i] = vec[i];
		for ( i = 1 ; i <= lg ; ++i ) {
			for ( strt = 0 ; strt < int(vec.size()) ; strt += ( 1 << (i+1) ) ) {
				md = strt + ( 1 << i );

				// prefijos
				if ( md < int(vec.size()) ) {
					acum = vec[md];
					table[i][md] = acum;
				}
				for ( j = md+1 ; j < strt+(1<<(i+1)) && j < int(vec.size()) ; ++j ) {
					acum = acum+ vec[j];
					table[i][j] = acum;
				}

				// sufijos
				if ( md-1 < int(vec.size()) ) {
					acum = vec[md-1];
					table[i][md-1] = acum;
				} else {
					acum = vec[vec.size()-1];
					table[i][vec.size()-1] = acum;
				}
				for ( j = min( md-2, vec.size()-2 ) ; j >= strt ; --j ) {
					acum = acum+ vec[j];
					table[i][j] = acum;
				}
			}
		}
	}
	Typ query( int l, int r ) {
		if ( l == r ) return table[0][l];
		int cut = ilog2( l^r );
		Typ res = table[cut][l]+ table[cut][r];
		return res;
	}
};