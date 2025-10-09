/*
Autor: Oscar Vargas Pabon

Esto funciona en O(nlgn) de memoria y preprocesamiento, O(lgn) de query

Notar que estoy asumiendo que 0 es la raiz (de lo contrario se puede danar 
	la acumulacion de bl "bl[nd][i+1]" en build_lca)

Fue testeado en: https://codeforces.com/problemset/problem/1843/F2
*/
class Data{
	public:
	int sm,mxpr,mnpr,mxsf,mnsf,mxsq,mnsq;
	Data(int xx=0){
		sm=xx;
		mxpr=mxsf=mxsq=max(0,xx);
		mnpr=mnsf=mnsq=min(0,xx);
	}
	Data operator +(const Data &o)const{
		Data neo;
		neo.sm=sm+o.sm;
		neo.mxpr=max(mxpr,sm+o.mxpr);
		neo.mnpr=min(mnpr,sm+o.mnpr);
		
		neo.mxsf=max(o.mxsf,o.sm+mxsf);
		neo.mnsf=min(o.mnsf,o.sm+mnsf);
		
		neo.mxsq=max({mxsq,o.mxsq,mxsf+o.mxpr});
		neo.mnsq=min({mnsq,o.mnsq,mnsf+o.mnpr});
		return neo;
	}
	void invert(){ swap(mnpr,mnsf); swap(mxpr,mxsf); }
};

const int lgi=20;
Data bl[template_limit][lgi];
int pi[template_limit][lgi], dpt[template_limit];

void build_lca(const vector<vector<int>> &t, const vector<Data> &arr, int nd=0,int p=0, int d=0 ){
	dpt[nd]=d;
	pi[nd][0]=p;
	bl[nd][0]=arr[nd];
	
	rep(i,0,lgi-1){
		int anc = pi[nd][i];
		pi[nd][i+1]=pi[anc][i];
		if(anc)bl[nd][i+1]=bl[nd][i] + bl[anc][i];
	}
	
	for(int e:t[nd])if(e!=p)build_lca(t,arr,e,nd,d+1);
}
Data query(int u,int v){
	auto jump=[&](Data &dt, int &x,int j){ dt=dt+bl[x][j]; x=pi[x][j]; };
	if(dpt[u]<dpt[v])swap(u,v);
	Data l,r;
	rep(i,lgi-1,-1)if( (dpt[u]-dpt[v])&(1<<i) ) jump(l,u,i);
	if(u==v)return l+bl[u][0];
	
	rep(i,lgi-1,-1)if( pi[u][i]!=pi[v][i] ){
		jump(l,u,i); jump(r,v,i);
	}
	jump(l,u,0); jump(r,v,0);
	l=l+bl[u][0]; r.invert();
	return l+r;
}
