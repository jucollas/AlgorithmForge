
#include <bits/stdc++.h>

typedef long long lint;

using namespace std;

#define rep(i,strt,end) for(int i = strt ; i !=int(end) ; (int(strt)<int(end))?++i:--i )
const lint mod =998244353;
/*
Autor: Oscar Vargas Pabon

and convolution tested in https://judge.yosupo.jp/problem/bitwise_and_convolution
or convolution currently untested

Note it wont work well for A*A

Taken from https://codeforces.com/blog/entry/119082
*/


template<typename T>
class btw_conv{
public:
	void resize(vector<T> &vec){
		int e=1;while(e<int(vec.size()))e*=2;
		vec.resize(e);
	}

	void trans_subset(vector<T> &vec,int sd){
		// I assume sd\in\{-1,1\}
		for(int e=1;e<int(vec.size());e*=2)rep(i,0,vec.size()){
			if(i&e)vec[i]+=vec[i^e]*sd,vec[i]%=mod;
		}
	}
	void or_conv(vector<T> &A,vector<T> &B){
		// the answer is shown in A; I assume |A|=|B|=2^x for some x
		trans_subset(A,1);trans_subset(B,1);
		rep(i,0,A.size())A[i]*=B[i];
		trans_subset(A,-1);
	}
	
	void trans_superset(vector<T> &vec,int sd){
		// I assume sd\in\{-1,1\}
		for(int e=1;e<int(vec.size());e*=2)rep(i,0,vec.size()){
			if(i&e)vec[i^e]+=vec[i]*sd,vec[i^e]%=mod;
		}
	}
	
	void and_conv(vector<T> &A,vector<T> &B){
		// the answer is shown in A; I assume |A|=|B|=2^x for some x
		trans_superset(A,1);trans_superset(B,1);
		rep(i,0,A.size())A[i]*=B[i];
		trans_superset(A,-1);
	}
};
btw_conv<lint> doall;


/*
Autor: Oscar Vargas Pabon

GCD_CONV Currently UNTESTED
LCM_CONV tested in https://codeforces.com/gym/105053/problem/G

Note it wont work well for A*A

Taken from https://codeforces.com/blog/entry/119082

Im assuming from my template:
#define rep(i,strt,end) for(int i = strt ; i !=int(end) ; (int(strt)<int(end))?++i:--i )
*/


template<typename T>
class div_conv{
	vector<int> prm;
public:
	div_conv(int limit){
		vector<bool> crb(limit,1);crb[0]=crb[1]=0;
		rep(i,0,limit)if(crb[i]){
			prm.push_back(i);
			if(i*1ll*i<lint(limit))for(int j=i*i;j<limit;j+=i){
				crb[j]=0;
			}
		}
	}

	void zeta_mult(vector<T> &vc){
		// vc'_x=\sum_{x|y}vc_y
		for(int p:prm)for(int i=(vc.size()-1)/p;i;--i){
			vc[i]+=vc[i*p]; vc[i]%=mod;
		}
	}
	void mobi_mult(vector<T> &vc){
		// vc_x=\sum_{x|y}vc'_y
		for(int p:prm)for(int i=1;i*p<int(vc.size());++i){
			vc[i]-=vc[i*p]; vc[i]%=mod;
		}
	}
	void gcd_conv(vector<T> &A,vector<T> &B){
		// I assume |A|=|B|
		// A'_x=\sum_{x=gcd(u,v)}A_uB_v
		zeta_mult(A);zeta_mult(B);
		rep(i,0,A.size())A[i]*=B[i];
		mobi_mult(A);
	}
	
	void zeta_div(vector<T> &vc){
		// vc_x=\sum_{y|x}vc'_y
		for(int p:prm)for(int i=1;i*p<int(vc.size());++i){
			vc[i*p]+=vc[i]; vc[i*p]%=mod;
		}
	}
	void mobi_div(vector<T> &vc){
		// vc'_x=\sum_{y|x}vc_y
		for(int p:prm)for(int i=(vc.size()-1)/p;i;--i){
			vc[i*p]-=vc[i]; vc[i*p] %=mod;
		}
	}
	
	void lcm_conv(vector<T> &A,vector<T> &B){
		// I assume |A|=|B|
		// A'_x=\sum_{x=lcm(u,v)}A_uB_v
		zeta_div(A);zeta_div(B);
		rep(i,0,A.size())A[i]*=B[i];
		mobi_div(A);
	}

};

div_conv<lint> dc(1e6+1);
// div_conv<lint> dc(10);

int main(){
	ios_base::sync_with_stdio(false);
    cin.tie(NULL);
	
	int n;cin>>n;
	/*vector<lint> a(1<<n),b(1<<n);
	rep(i,0,1<<n)cin>>a[i];
	rep(i,0,1<<n)cin>>b[i];
	
	doall.and_conv(a,b);*/
	vector<lint>a(n+1),b(n+1);
	rep(i,1,n+1)cin>>a[i];
	rep(i,1,n+1)cin>>b[i];
	
	// dc.zeta_mult(a);
	// repcout << ac << ' ';cout << endl;
	// dc.mobi_mult(a);
	// for(lint ac:a)cout << ac << ' ';cout << endl;
	
	// dc.gcd_conv(a,b);
	dc.lcm_conv(a,b);
	rep(i,1,n+1)cout << ((a[i]%mod)+mod)%mod << ' ';
	cout <<'\n';
	
	return 0;
}