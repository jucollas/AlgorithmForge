/*
Autor: Oscar Vargas Pabon
Construyo en O(nlgn) con memoria O(nlg n)

Requiero conmutatividad
query en O(1) pero ademas requiero idempotencia
binlift en O(lg n)
*/
// int ilog2( int num ) { return (sizeof(int) * 8) - __builtin_clz( num ) - 1; }
// asumo lo anterior de mi template
typedef int Typ;
class SparseTable {
public:
	vector<vector<Typ>> table;
	SparseTable() = default;
	SparseTable( const vector<Typ> &vec ) {
		int lg = ilog2( vec.size() ), i, j;
		table = vector<vector<Typ>> ( lg+1, vector<Typ> ( vec.size() )  );
		for ( i = 0 ; i < int(vec.size()) ; ++i ) table[0][i] = vec[i];
		for ( i = 1 ; i <= lg ; ++i ) {
			for ( j = 0 ; j + ( 1 << i ) <= int(vec.size()) ; ++j ) {
				table[i][j] = table[i-1][j] + table[i-1][j + (1<<(i-1)) ];
			}
		}
	}
	Typ query( int l, int r ) {
		// requiere idempotencia
		int lg = ilog2( r-l+1 );
		Typ res =table[lg][l] + table[lg][r- (1<<(lg)) +1 ];
		return res;
	}
	Typ binlift( int l, int r ) {
		Typ res;
		for ( int i = table.size()-1 ; i >= 0 && l <= r ; --i ) {
			if ( l + ( 1 << i )-1 <= r ) {
				res = res + table[i][l];
				l += ( 1 << i );
			}
		}
		return res;
	}
};