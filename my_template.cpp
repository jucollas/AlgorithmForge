/*
Autor: Oscar Vargas Pabon
Fecha: 

*/

#include <bits/stdc++.h>

typedef long long lint;

using namespace std;

#define rep(i,strt,end) for(int i = strt ; i !=int(end) ; (int(strt)<int(end))?++i:--i )
#define rall(vec) vec.rbegin(), vec.rend()
#define all(vec) vec.begin(), vec.end()
#define pb push_back
#define ppb pop_back
#define pf push_front
#define ppf pop_front

mt19937_64 random_64( chrono::steady_clock::now().time_since_epoch().count() );
int ilog2( int num ) { return 8*sizeof(int) - __builtin_clz( num ) - 1; }
lint mpow(lint x,lint e,lint m){lint res=1ll;while(e){if(e&1ll)res=(res*x)%m;e>>=1;x=(x*x)%m;}return res;}

const int template_limit = 1e6;
int a[template_limit], b[template_limit];

void solve() {
	
}

int32_t main(){
	ios_base::sync_with_stdio(false);
    cin.tie(NULL);
	cout << setprecision(12) << fixed;

    int t = 2;
    cin >> t; ++t;
    while ( --t ) {
		solve();
    }
	return 0;
}

