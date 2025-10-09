/*
Probando mis implementaciones para la ICPC


*/

#include <bits/stdc++.h>

using namespace std;
#define pb push_back

/* Insertar aqui */

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

/* El resto _-_-_- */

vector<int> slowVersion( const string &str ) {
	vector<int> res( str.size(), 0 );
	for ( int i = 0 ; i < int(str.size()) ; ++i ) {
		for ( int ran = 0 ; ran <= i ; ++ran ) {
			bool tmp = 1;
			for ( int j = 0 ; j < ran && tmp ; ++j ) tmp = str[j]==str[i-ran+1+j];
			if ( tmp ) res[i]=ran;	
		}
	}
	return res;
}

bool is_equal( const vector<int> &a, const vector<int> &b ) {
	bool res = a.size()==b.size();
	for ( int i = 0 ; i < int(a.size()) && res ; ++i ) res = a[i] == b[i];
	return res;
}

int random( int tope ) {
    int i = abs( rand() % tope );
    return i;
}
void printArr( const vector<int> &arr, const string &str ) {
	cout << "[" << str << "]";
	for ( int i = 0 ; i < int(arr.size()) ; ++i ) cout << ' ' <<arr[i];
	cout << endl;
}

int main() {
	const int sz=1000000, tope = 1<<7;
	while ( 1 ) {
		string cad;
		for ( int i = 0 ; i < sz ; ++i ) cad.pb( random( tope )+'a' );
		vector<int> v1 = slowVersion( cad );
		vector<int> v2 = prefixFunction(cad );
		if ( !is_equal( v1, v2 ) ) {
			cout << cad << endl;
			printArr( v1, "slow" );
			printArr(v2,"impl" );
		}
		assert ( is_equal(v1,v2) );
	}
	return 0;
}