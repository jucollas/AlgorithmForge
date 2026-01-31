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
void raw_debug(istream_iterator<string> it) {cerr<<endl;assert(it==it);}
template<typename T, typename... Args>
void raw_debug(istream_iterator<string> it, T a, Args... args) { cerr <<"<"<< *it << "->" << a << "> "; raw_debug(++it, args...); }
#define idebug(v) {cerr<<'['<<#v<<']';for(const auto &el:v)cerr << ' ' << el; cerr << endl;}
#define adebug(ar,n) {cerr<<'['<<#ar<<']';for(int i=0;i<n;++i)cerr << ' ' << ar[i]; cerr << endl;}

#define rep(i,strt,end) for(int i = strt ; i !=int(end) ; (int(strt)<int(end))?++i:--i )
#define rall(vec) vec.rbegin(), vec.rend()
#define all(vec) vec.begin(), vec.end()
#define sz(vec) int(vec.size())
#define eb emplace_back
#define pb push_back
#define pob pop_back
#define pf push_front
#define pof pop_front

lint seed=chrono::steady_clock::now().time_since_epoch().count();
// lint seed=1769878272471046000;
mt19937_64 rng_64( seed );
constexpr int ilog2( int num ) { return 8*sizeof(int) - __builtin_clz( num ) - 1; }
template<typename tpow> constexpr tpow mpow(tpow x,lint e,tpow m){tpow res=1;while(e){if(e&1ll)res=(res*1ll*x)%m;e>>=1;x=(x*1ll*x)%m;}return res;}

const int mod=1e9+7;
#include"..\modulo_int.cpp"
#include"..\set_power_series.cpp"

template<typename T> vector<T> brute_or_conv(const vector<T> &A,const vector<T>&B){
	const int n=A.size();
	vector<T> C(A.size(),0);rep(i,0,n)rep(j,0,n)C[i|j]+=A[i]*B[j];
	return C;
}
template<typename T> vector<T> brute_and_conv(const vector<T> &A,const vector<T>&B){
	const int n=A.size();
	vector<T> C(A.size(),0);rep(i,0,n)rep(j,0,n)C[i&j]+=A[i]*B[j];
	return C;
}
template<typename T> vector<T> brute_sub_conv(const vector<T> &A,const vector<T>&B){
	const int n=A.size();
	vector<T> C(A.size(),0);rep(i,0,n)rep(j,0,n)C[i|j]+=(A[i]*B[j])*T((i&j)==0);
	return C;
}
template<typename T> vector<T> brute_set_exp(const vector<T> &A){
	assert(A[0]==0);
	vector<T> xd(A.size(),0),Ai(A.size(),0);Ai[0]=1;
	mint fact=1,cnt=1;
	for(int e=1;e<=int(A.size());e<<=1){
		rep(i,0,A.size())xd[i]+=Ai[i]/fact;
		Ai=subset_conv<T>(Ai,A);
		fact*=cnt;cnt+=1;
	}return xd;
}
void solve() {
	const int lgi=10,iter=1e4;
	rep(tm,0,iter){
		debug(tm,iter);
		vector<mint> A(1<<lgi),B=A,BB;
		rep(i,0,1<<lgi)A[i]=rng_64(),B[i]=rng_64();
		while(A[0]==0)A[0]=rng_64();
		
		vector<mint> b_or,r_or,b_and,r_and,b_sub,r_sub,b_exp,r_exp,r_log;
		
		b_or=brute_or_conv<mint>(A,B);
		r_or=or_conv<mint>(A,B);
		
		bool or_ok=b_or==r_or;
		// debug("or",or_ok);
		
		b_and=brute_and_conv<mint>(A,B);
		r_and=and_conv<mint>(A,B);
		
		bool and_ok=b_and==r_and;
		// debug("and",and_ok);
		
		b_sub=brute_sub_conv<mint>(A,B);
		r_sub=subset_conv<mint>(A,B);
		// b_sub=r_sub;
		
		bool sub_ok=b_sub==r_sub;
		// debug("sub",sub_ok);
		
		vector<mint> r_isub=subset_iconv<mint>(A,r_sub);
		bool isub_ok=r_isub==B;
		
		mint tmp=A[0];A[0]=0;
		b_exp=brute_set_exp<mint>(A);
		r_exp=set_exp<mint>(A);
		bool exp_ok=b_exp==r_exp;
		
		r_log=A;
		r_log=set_log<mint>(r_exp);
		bool log_ok=A==r_log;
		
		A[0]=tmp;
		
		
		
		
		if(!or_ok||!and_ok||!sub_ok||!isub_ok||!exp_ok||!log_ok){
			debug(or_ok,and_ok,sub_ok,isub_ok,exp_ok,log_ok);
			idebug(A);idebug(B);
			idebug(r_isub);
			
			debug(seed);
			
			
			// idebug(b_sub);idebug(r_sub);
			A[0]=0;idebug(A);idebug(r_exp);idebug(r_log);
			return;
		}
	}debug("ok");
}

int32_t main(){
	ios_base::sync_with_stdio(false);
    cin.tie(NULL);
	cout << setprecision(12) << fixed;

    int t = 2;
    // cin >> t; ++t;
    while ( --t ) {
		solve();
    }
	return 0;
}

