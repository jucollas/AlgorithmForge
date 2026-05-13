/*
Author:Oscar Vargas Pabon
Computing: Given $A[0..n)$ and $B[0..m)$ compute
	$B[0..n+m-1)$ where $B[k]=max_{i+j=k}\{A[i]+B[j]\}$
The concave version is basically minkowski sum of convex polygons 
	but specialized for this specific problem.
Concave version tested in: https://atcoder.jp/contests/ndpc/tasks/ndpc2026_o

Ill eventually keep looking 
*/
template<typename tint>
void assert_concavity(const vector<tint>&A){
	if(A.size()<=1)return; // O(n)
	tint prv=A[1]-A[0];rep(i,2,A.size()){
		const tint delta=A[i]-A[i-1];
		if(delta>prv)assert(string("Esto no es")==string("concavo"));
		if(delta<prv)assert(string("Esto no es")==string("convexo"));
		prv=delta;
	} }
template<typename tint>
vector<tint> maxplus_convolution_concave(const vector<tint>&A,const vector<tint>&B){
	const int n=A.size(),m=B.size();//O(n+m)
	vector<tint> res(n+m-1);
	int aind=0,bind=0;while(aind<n&&bind<m){
		res[aind+bind]=A[aind]+B[bind];
		if(aind+1>=n)++bind;
		else if(bind+1>=m)++aind;
		else if(A[aind+1]-A[aind]>=B[bind+1]-B[bind])++aind;
		else ++bind;
	} return res; }
template<typename tint,tint inf=1e9>
vector<tint> max_plus_convolution_naive(){//O(n*m)
	vector<tint> brute(n+m-1,-inf);rep(i,0,n)rep(j,0,m)
		brute[i+j]=max<tint>(brute[i+j],A[i]+B[j]);
	return brute; }