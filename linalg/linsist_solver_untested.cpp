/*
Autor:Oscar Vargas Pabon
*/

typedef int tmat;
typedef vector<vector<tmat>> mat;
vector<int> plu_decomp( const mat&m, mat &l, mat&r ){
	assert(!m.empty()); assert(m.size()==m[0].size());
	int n=m.size();
	vector<int> pi(n);//permutacion
	l=mat(n,vector<tmat>(n,0));r=m;
	rep(i,0,n) l[i][i]=1;//inicio l
	
	rep(i,0,n) {
		
		int nr=i; // hallo el pivote
		rep(j,i+1,n) if(abs(u[nr][i])<abs(u[j][i]))nr=j;
		swap(pi[i],pi[nr]);
		rep(j,0,n)swap(m[i][j],m[nr][j]);
		
		rep(j,i+1,n){ // eliminacion gaussiana
			l[j][i] = u[j][i]/u[i][i];
			rep(k,i,n) u[j][k] -= l[j][i]*u[i][k];
		}
	}
	return pi;
}
vector<tmat> solve_diag( const vector<int> &pi, const mat&m, const vector<tmat>&y, bool l){
	// m*x=y ; output -> x ; l indica si es lower(l==1) o upper diagonal(l==0)
	int n = y.size();
	vector<tmat> x(n,0);
	if ( l ) {
		rep(i,0,n) {
			rep(j,0,i) x[pi[i]]+= x[pi[j]]*m[i][j];
			x[pi[i]] = (y[pi[i]]-x[pi[i]])/m[i][i];
		}
	} else {
		rep(i,n-1,-1){
			rep(j,i+1,n) x[pi[i]]+= x[pi[j]]*m[i][j];
			x[pi[i]] = (y[pi[i]]-x[pi[i]])/m[i][i];
		}
	}
	return x;
}
vector<tmat> plu_solve(const mat&m,const vector<tmat> &y){
	// soluciona m*x=y ; retorna x
	mat l,r;
	vector<int> pi=plu_decomp(m,l,r);
	vector<tmat> x = solve_diag(pi,l,solve_diag(pi,r,y,0),1);
	return x;
}
vector<tmat> simple_iteration(const mat&m,const vector<tmat> &y,int iter=10){
	// hallo Ax=y ; luego itero sobre Ax'=y-Ax, haciendo luego x=x+x'
	int n=m.size();
	vector<tmat> x = plu_solve(m,y),xe,ye;
	rep(i,0,iter){
		rep(i,0,n) {
			ye[i]=y[i];
			rep(j,0,n)ye[i]-=m[i][j]*x[j];
		}
		xe = plu_solve(m,ye);
		rep(i,0,n) x[i]+=xe[i];
	}
	return x;
}