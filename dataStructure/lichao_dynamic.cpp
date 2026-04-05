// Author: Oscar Vargas Pabon. Motivated by Judival
// tested in https://codeforces.com/problemset/problem/932/F
// requires Point from ..\geometry\geotemplate.cpp supporting dot product (*)
// and '=='.
// It works by assigning an id to each tree. It supports update, alias add line;
// querys and merge trees. The tree has size O(n) for n:linear functions added
// and query/update time O(lgn), merge time O(nlgm) for m being the tree being merged to
template<typename Tpt,lint tl=lint(-1e9),lint tr=lint(1e9),Tpt inf=Tpt(1e18)>
struct LiChao_dyn{
	static Tpt eval(const Point<Tpt>&ln,Tpt x){return Point<Tpt>(x,1)*ln;}
	struct Node{
		Point<Tpt> ln; array<int,2> chld;
		void reset(){ln={0,inf};chld={-1,-1};}
		Node(Point<Tpt> aln={0,inf}){reset();ln=aln;}
	};
	vector<Node> pool;int pool_sz; void reset(){pool_sz=0;}
	constexpr LiChao_dyn(int preres=0){reset();pool.reserve(preres);}
	int new_tree(){if(int(pool.size())==pool_sz)pool.push_back(Node());
		 pool[pool_sz].reset(); int rs=pool_sz++; return rs;}
 
	int update(int nd, Point<Tpt>ln, lint l=tl,lint r=tr){
		if(nd==-1) nd=new_tree();
		const lint m=(l+r)/2ll; Point<Tpt>&aln=pool[nd].ln;
		if(eval(ln,m)<eval(aln,m))swap(ln,aln);
		if(ln==Point<Tpt>(0,inf)||l>=r)return nd;
		if     (eval(ln,l)<eval(aln,l))pool[nd].chld[0]=update(pool[nd].chld[0],ln, l ,m);
		else if(eval(ln,r)<eval(aln,r))pool[nd].chld[1]=update(pool[nd].chld[1],ln,m+1,r);
		return nd;
	} Tpt query(int nd,Tpt x, lint l=tl,lint r=tr){
		if(nd==-1)return inf;
		const lint m=(l+r)/2ll; Tpt res=eval(pool[nd].ln,x);
		if(x<=m) res=min<Tpt>(res,query(pool[nd].chld[0],x, l ,m));
		else     res=min<Tpt>(res,query(pool[nd].chld[1],x,m+1,r));
		return res;
	} int merge(int to,int from,lint l=tl,lint r=tr){
		if(to==-1)return from;
		else if(from==-1)return to;
		to=update(to,pool[from].ln,l,r); const lint m=(l+r)/2ll;
		const array<lint,2> bndl={l,m+1},bndr={m,r}; for(int i=0;i<2;++i)
			pool[to].chld[i]=merge(pool[to].chld[i],pool[from].chld[i],bndl[i],bndr[i]);
		return to;
	}
}; LiChao_dyn<lint,lint(-1e5-33),lint(1e5+33)> cht(1e5);