/*
Autor: Oscar Vargas Pabon
Fecha: 25/04/2024

rmq using sparse table in <O(n), O(1)>
https://codeforces.com/blg/entry/78931

Nota: cambiar 'myMin' por lo mismo mediante operador ternario puede mejorar un poco
		el tiempo de ejecucion segun el blg
*/
//int ilog2( int num ) { return sizeof(int)*8 -__builtin_clz(num)-1; }
// asumo lo anterior de mi template
typedef int Typ;
class rmqSparseTable {
	private:
	static const int b = 30;
	public:
	vector<Typ> table, mask, arr;
	
	int n;
	
	static Typ myMin( Typ n1, Typ n2 ) {
		Typ res = ( n1 < n2 ) ? n1 : n2;
		return res;
	}
	
	Typ queryBlock( int r, int sz = b ) {
		int lg = ilog2( mask[r] & ((1<<sz)-1) );
		return arr[r-lg];
	}
	
	rmqSparseTable( const vector<Typ> &arrIn ) {
		arr = arrIn;
		n = arr.size();
		table.resize( n );
		mask.resize( n );
		
		int actm = 0;
		for ( int i = 0 ; i < n ; ++i ) {
			// lo limito a b bits
			actm = ( actm << 1 ) & ( ( 1 << b )-1 );
			// quitamos de la cola monotonica implicita en actm
			while ( actm > 0 && arr[i] <= arr[i-ilog2(actm&(-actm))] ) {
				actm ^= actm&(-actm);
			}
			actm |= 1;
			mask[i] = actm;
		}
		//llenamos la sparse table
		for ( int i = 0 ; i < n/b ; ++i ) table[i] = queryBlock( (i+1)*b -1 );
		for ( int j = 1 ; (1<<j) <= n/b ; ++j ) {
			for ( int i = 0 ; (1<<j)+i <= n/b ; ++i ) {
				table[(n/b)*j + i] = myMin( table[(n/b)*(j-1) + i], table[(n/b)*(j-1)+ i + (1<<(j-1))] );
			}
		}
	}
	
	Typ query( int l, int r ) {
		Typ res;
		if ( r-l+1 <= b ) {
			res = queryBlock( r, r-l+1 );
		} else {
			res = myMin( queryBlock(l+b-1), queryBlock(r) );
			
			// se añade el +1 porque no queremos que encuentre el menor en el bloque en donde esta l porque puede salirse del rango
			int neol = l/b+1, neor = r/b-1;
			if ( neol <= neor ) {
				int lg = ilog2( neor-neol+1 );
				res = myMin( res, myMin( table[(n/b)*lg+neol], table[(n/b)*lg+neor-(1<<lg)+1] ) );
			}
		}
		return res;
	}
};