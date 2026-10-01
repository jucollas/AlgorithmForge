/* Author: Oscar Vargas Pabon
Some random tree-tricks
*/  ////// this trick was first used in
// https://codeforces.com/gym/102201/problem/F
int tin[max_n],tout[max_n]; // taking in-out times (my version)
int get_time(const vector<vector<int>>&t,int nd=0,int p=0,int tm=0){
	tin[nd]=tm++;for(int e:t[nd])if(e!=p)tm=get_time(t,e,nd,tm);
	return tout[nd]=tm;
} // when having in and out times like this, one can query path sums
int a[max_n],acm[max_n]; // in the following way
void built(int n){ // the indexing is a bit fragile, so take that into
	for(int i=0;i<=n;++i)acm[i+1]=0; // account. 
	for(int i=0;i<n;++i)acm[tin[i]+1]+=a[i],acm[tout[i]+1]-=a[i];
	for(int i=0;i<n;++i)acm[i+1]+=acm[i];
} inline int root_sum(int nd){ return acm[tin[nd]+1]; }
inline int path_sum(int u,int v){
	return root_sum(u)+root_sum(v)-2*root_sum(lca(u,v))+a[lca(u,v)];
} // in Some cases, it can be better to do this via dfs and mantaining
// some shared memory. The idea is to keep at all times the answers
// from the root to the actual node
