#include<bits/stdc++.h>
using namespace std;



// salio de https://codeforces.com/blog/entry/15643
#define debug(args...) { string _s = #args; replace(_s.begin(), _s.end(), ',', ' '); stringstream _ss(_s); istream_iterator<string> _it(_ss); raw_debug(_it, args);}
void raw_debug(istream_iterator<string> it) {cerr<<endl;}
template<typename T, typename... Args>
void raw_debug(istream_iterator<string> it, T a, Args... args) { cerr <<"<"<< *it << "->" << a << "> "; raw_debug(++it, args...); }

// #define t1(a) cout << #a << endl;
// #define t1(a,args...) t1(args)
//#define hola(act,args...) {cout << act << ' ';hola(args);}
#define t2(args...) string _hh=#args;cout<<_hh<<endl;
//#define rep(i,strt,end) for(__typeof(strt) i = strt ; i !=end ; (strt<end)?++i:--i )
#define rep(i,strt,end) for(__typeof(strt) i = strt ; i !=(__typeof(strt))end ; (strt<(__typeof(strt))end)?++i:--i )
#define all(str) str.begin(),str.end()


// #define dball2(a,b) cout << *a << endl;
// #define dball1(v) dball2(all(v),all(v))

#define xx (4 + 1)
#define yy (2 * xx)

#define test(aaa)cout << #aaa << endl;
#define f(a,b) cout << a##b<<endl;
int main(){
	rep(i,0,10)cout << i << ' '; cout << endl;
	rep(i,10,0)cout << i << ' '; cout << endl;
	int x,y;
	pair<int,int> a(1,2);
	tie(x,y)=a;
	int z;
	tie(ignore,z)=a;
	cout << x << " _ " << y << " _ " << z << endl;
	cout << "pasamos aca" << endl;
	t2(x,y,z);
	// hola(x,y,z);
	// t1(x,y,z);
	debug(x,y,z);
	string hola="degenere";
	debug(hola,"adios");
	debug(a.first,a.second);
	vector<int> b={8,1,2,3,4,5,6};
	list<int> c={1,2,3,4,5,6};
	set<int> d = {1,2,3,4,5,6};
	sort(all(b));
	test(all(b));
	rep(it,b.begin(),b.end()) cout << *it << ' ';
	// rep(it,all(b))cout << *it << ' ';
	// rep(it,all(c))cout << *it << ' ';
	// rep(it,all(d))cout << *it << ' ';
	// f(1,"h");
	// dball1(b);
	cout << endl;
	cout << yy << "  __ _ " << endl;
	return 0;
}
