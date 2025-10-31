/*
Autor: Oscar Vargas Pabon

mat_mult tested in UVA 13298 (A Fibonacci Family Formula)
*/

typedef int tmat;
typedef vector<vector<tmat>> mat;
mat mat_mult(const mat &m1,const mat &m2){
	// O(n*m*h)
	assert(!m2.empty()); assert(!m1.empty()); assert(m2.size()==m1[0].size());
	int n = m1.size(),m=m2.size(),h=m2[0].size();
	mat res(n,vector<tmat>(h,0));
	// rep(i,0,n) rep(j,0,h) rep(k,0,m) res[i][j] += m1[i][k]*m2[k][j];
	rep(i,0,n) rep(j,0,h) rep(k,0,m) res[i][j] = (res[i][j]+ m1[i][k]*1ll*m2[k][j])%mod;
	return res;
}
