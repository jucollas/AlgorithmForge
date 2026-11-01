/* Author: claude
Test de manacher
Build: g++ -std=c++17 -O1 -g -D_GLIBCXX_ASSERTIONS -fsanitize=address,undefined manacherTest.cpp -o /tmp/mt && /tmp/mt
Para cada posicion i de c'="&#a#b#...#%" calcula por fuerza bruta R[i], el maximo x tal que
c'[i-x..i+x] es palindromo (sin tocar los centinelas). R[i] es justamente la longitud del
palindromo mas largo de la cadena original centrado ahi.
Se compara contra lo que promete el comentario: res[i]==R[i].
*/
#include <bits/stdc++.h>
using namespace std;

#include "../manacher_notMine.cpp"

string transformed( const string &cad ) {
	string c2="&#"; for(char c:cad){c2.push_back(c);c2.push_back('#');} c2.push_back('%');
	return c2;
}
vector<int> slowVersion( const string &cad ) {
	string c2=transformed(cad); int n=c2.size(); vector<int> R(n-1,-1);
	for ( int i = 1 ; i < n-1 ; ++i ) {
		R[i]=0;
		while ( c2[i-R[i]-1]==c2[i+R[i]+1] ) ++R[i];
	}for(int i=0;i+2<n-1;++i)R[i]=R[i+2];
	R.pop_back();R.pop_back();if(!R.empty())R.pop_back();
	return R;
}

mt19937 rng(1234);
string randomString( int n, int sigma ) {
	string s(n,'a'); for(char &c:s) c='a'+rng()%sigma; return s;
}
bool check( const string &s ) {
	// printf("what the fuck\n");
	// cout << s << " -> " << s.size() << endl;

	vector<int> want=slowVersion(s),got=manacher(s);
	// cerr << want.size()<< " _ " << got.size() << endl;

	// cerr << "got";for(int i=0;i<got.size();++i)cerr << ' ' << got[i];cerr << endl;
	// cerr << "wnt";for(int i=0;i<want.size();++i)cerr << ' ' << want[i];cerr << endl;
	
	// printf("nemompaf\n");
	if ( got.size()!=want.size() ) { 
		printf("FAIL \"%s\": size %zu, expected %zu\n",s.c_str(),got.size(),want.size());
		// exit(0);
		return false;
	}
	for ( int i = 1 ; i < int(want.size()) ; ++i ) if ( got[i]!=want[i] ) {
		printf("FAIL \"%s\" (c'=\"%s\"): res[%d]=%d, expected %d\n",s.c_str(),transformed(s).c_str(),i,got[i],want[i]);
		printf("  got: "); for(int i=1;i<int(got.size());++i)printf("%d ",got[i]);
		printf("\n want: "); for(int i=1;i<int(want.size());++i)printf("%d ",want[i]); printf("\n");
		// exit(0);
		return false;
	}
	return true;
}

int main() {
	// check("aa");
	int fails=0, total=0;
	for ( string s : {"","a","aa","ab","aba","abba","aaaa","abacaba","abcba","banana"} ) fails+=!check(s), ++total;
	// todas las cadenas sobre {a,b} hasta longitud 12
	// printf("??????");
	for ( int n = 1 ; n <= 12 && fails < 5 ; ++n ) for ( int mask = 0 ; mask < (1<<n) && fails < 5 ; ++mask ) {
		string s(n,'a'); for(int j=0;j<n;++j) if(mask>>j&1) s[j]='b';
		fails+=!check(s), ++total;
	}
	// printf("lleha aca \n");
	for ( int it = 0 ; it < 20000 && fails < 5 ; ++it ) fails+=!check(randomString(rng()%200+1,rng()%4+1)), ++total;
	for ( int it = 0 ; it < 5 && fails < 5 ; ++it ) fails+=!check(randomString(100000,2)), ++total; // n grande (bruta es O(n^2) en el peor caso, pero ok con sigma=2 aleatorio)
	if ( fails ) printf("%d failing case(s) (stopped early) Manacher\n",fails);
	else printf("OK (%d strings) Manacher\n",total);
	return fails!=0;
}
