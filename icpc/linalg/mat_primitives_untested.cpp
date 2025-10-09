/*
Autor: Oscar Vargas Pabon
*/

typedef int tmat;
typedef vector<vector<tmat>> mat;
mat mult(const mat &m1,const mat &m2){
	// O(n*m*h)
	assert(!m2.empty()); assert(!m1.empty()); assert(m2.size()==m1[0].size());
	int n = m1.size(),m=m2.size(),h=m2[0].size();
	mat res(n,vector<tmat>(h,0));
	rep(i,0,n) rep(j,0,h) rep(k,0,m) res[i][j] += m1[i][k]*m2[k][j];
	return res;
}
