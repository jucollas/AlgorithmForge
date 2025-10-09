/*
Autor: Oscar Vargas Pabon

No olvidar el guardar las nuevas versiones en rts[rt_n++]
Esta version no tiene LAZY PROPAGATION

Probado en https://www.spoj.com/problems/DQUERY/
*/

class Data {
	public:
	int sum;
	Data(int s=0) : sum(s) {};
	Data operator + ( const Data &o ) const {
		return Data(sum+o.sum);
	}
};

class Node {
	public:
	Node *l,*r;
	Data dat;
	Node(){ l=r=NULL; }
};

int t_sz, rt_n;
Node *rts[template_limit]; // roots;
list<Node*> to_erase; // guarda todo lo que se debe borrar;
Node* build( const vector<Data> &arr, int l=0,int r=-1 ) {
	if ( r==-1 ) t_sz=r=arr.size()-1;
	Node *nd = new Node(); to_erase.push_back(nd);
	if ( l>=r ) nd->dat = arr[l];
	else {
		int m = (l+r)/2;
		nd->l = build(arr,l,m);
		nd->r = build(arr,m+1,r);
		nd->dat = nd->l->dat + nd->r->dat;
	}
	return nd;
}
Node *update( int x, const Data &vl, Node *nd, int l=0, int r=t_sz ){
	Node *neo = new Node(); to_erase.push_back(neo);
	if ( l >= r ) neo->dat = vl;
	else {
		int m =(l+r)/2;
		if ( x <= m ) neo->l=update(x,vl,nd->l,l,m),neo->r=nd->r;
		else neo->r=update(x,vl,nd->r,m+1,r),neo->l=nd->l;
		neo->dat = neo->l->dat + neo->r->dat;
	}
	return neo;
}
Data query( int ql, int qr, Node *nd, int l=0, int r=t_sz ) {
	Data res;
	if ( ql<=l && r <= qr ) res = nd->dat;
	else if ( !(r<ql||qr<l) ) {
		int m = (l+r)/2;
		res = query(ql,qr,nd->l,l,m)+query(ql,qr,nd->r,m+1,r);
	}
	return res;
}

void free_pseg(){
	for ( Node *nd : to_erase) delete nd;
	rt_n = 0; to_erase={};
}