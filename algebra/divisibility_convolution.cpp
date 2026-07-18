/* Author: Oscar Vargas Pabon
GCD_CONV tested in https://judge.yosupo.jp/problem/gcd_convolution
LCM_CONV tested in https://codeforces.com/gym/105053/problem/G
				   https://judge.yosupo.jp/problem/lcm_convolution
Note it wont work well for A*A

Taken from https://codeforces.com/blog/entry/119082

Im assuming from my template:
#define rep(i,strt,end) for(int i = strt ; i !=int(end) ; (int(strt)<int(end))?++i:--i )
*/
template<typename T> struct div_conv{
	vector<int> prm; constexpr div_conv(int limit){
		vector<bool> crb(limit,1);
		rep(i,2,limit){ if(crb[i])prm.push_back(i);
			for(int p:prm){if(p*1ll*i>=limit)break;
				crb[i*p]=1;if(i%p==0)break; } }
	} void zeta_mult(vector<T> &vc)const{
		// $vc'_x=\sum_{x|y}vc_y$
		const int n=(vc.size()-1);
		for(int p:prm){ if(p>n)break;
			for(int i=n/p;i;--i) vc[i]+=vc[i*p]; }
	} void mobi_mult(vector<T> &vc)const{
		// $vc_x=\sum_{x|y}vc'_y$
		const int n=(vc.size()-1);
		for(int p:prm){ if(p>n)break;
			for(int i=1;i*p<=n;++i)vc[i]-=vc[i*p]; }
	} void gcd_conv(vector<T> &A,vector<T> &B)const{
		// I assume |A|=|B|
		// $A'_x=\sum_{x=gcd(u,v)}A_uB_v$
		zeta_mult(A);zeta_mult(B);
		rep(i,0,A.size())A[i]*=B[i];
		mobi_mult(A);
	} void zeta_div(vector<T> &vc)const{
		// $vc_x=\sum_{y|x}vc'_y$
		const int n=(vc.size()-1);
		for(int p:prm){ if(p>n)break;
			for(int i=1;i*p<=n;++i)vc[i*p]+=vc[i]; }
	} void mobi_div(vector<T> &vc)const{
		// $vc'_x=\sum_{y|x}vc_y$
		const int n=(vc.size()-1);
		for(int p:prm){ if(p>n)break;
			for(int i=n/p;i;--i) vc[i*p]-=vc[i]; }
	} void lcm_conv(vector<T> &A,vector<T> &B)const{
		// I assume |A|=|B|
		// $A'_x=\sum_{x=lcm(u,v)}A_uB_v$
		zeta_div(A);zeta_div(B);
		rep(i,0,A.size())A[i]*=B[i];
		mobi_div(A); }
}; const div_conv<mint> dc(1e6);