/*
Probando mis implementaciones para la ICPC


*/
#include<cassert>
#include <bits/stdc++.h>

using namespace std;
#define pb push_back

#include "../prefixFunction.cpp"

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

int main() { int t=100;
	const int sz=1e3, tope = 1<<7;
	while (t--) {
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
	}cerr << "termino bien " << endl;
	return 0;
}