/*
Probando mis implementaciones para la ICPC


*/

#include <bits/stdc++.h>

using namespace std;
#define pb push_back

/* Insertar aqui */

vector<int> zFunction( const string &cad ) {
/* Computa la funcion z en O(n). Esta es z[i] -> maximo prefijo comun de cad y cad[i...] */
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

/* El resto _-_-_- */

vector<int> slowVersion( const string &str ) {
	vector<int> res( str.size(), 0 );
	for ( int i = 1 ; i < int(str.size()) ; ++i ) {
		while ( i+res[i] < int(str.size()) && str[res[i]] == str[i+res[i]] ) ++res[i];
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
	const int sz=100000, tope = 1<<7;
	while ( true ) {
		string cad;
		for ( int i = 0 ; i < sz ; ++i ) cad.pb( random( tope )+'a' );
		vector<int> v1 = slowVersion( cad );
		vector<int> v2 = zFunction(cad );
		if ( !is_equal( v1, v2 ) ) {
			cout << cad << endl;
			printArr( v1, "slow" );
			printArr(v2,"impl" );
		}
		assert ( is_equal(v1,v2) );
	}
	return 0;
}