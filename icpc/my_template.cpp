/*
 ________
|    ___ |
|  ,',.(`|
| :  `'  |
| :) _  (|
|  `:_)_,|
|________|

Autor: Oscar Vargas Pabon
Fecha: 

*/

#include <bits/stdc++.h>

typedef long long lint;

using namespace std;

#define debug(args...) { string _s = #args; replace(_s.begin(), _s.end(), ',', ' '); stringstream _ss(_s); istream_iterator<string> _it(_ss); raw_debug(_it, args);}
void raw_debug(istream_iterator<string> it) {cerr<<endl;}
template<typename T, typename... Args>
void raw_debug(istream_iterator<string> it, T a, Args... args) { cerr <<"<"<< *it << "->" << a << "> "; raw_debug(++it, args...); }
#define idebug(v) {cout<<'['<<#v<<']';for(const auto &el:v)cout << ' ' << el; cout << endl;}
#define adebug(ar,n) {cout<<'['<<#ar<<']';for(int i=0;i<n;++i)cout << ' ' << ar[i]; cout << endl;}

#define rep(i,strt,end) for(int i = strt ; i !=int(end) ; (int(strt)<int(end))?++i:--i )
#define rall(vec) vec.rbegin(), vec.rend()
#define all(vec) vec.begin(), vec.end()
#define pb push_back
#define pob pop_back
#define pf push_front
#define pof pop_front

mt19937_64 rng_64( chrono::steady_clock::now().time_since_epoch().count() );
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

