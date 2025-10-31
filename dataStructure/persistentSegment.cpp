/*
Autor: Oscar Vargas Pabon

No olvidar el guardar las nuevas versiones en rts[rt_n++]
Esta version no tiene LAZY PROPAGATION

Probado en https://www.spoj.com/problems/DQUERY/
           https://codeforces.com/contest/2111/problem/G
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
vector<Node*> to_erase; // guarda todo lo que se debe borrar;
Node* create_node(){Node *nd=new Node();to_erase.push_back(nd);return nd;}
Node* build( const vector<Data> &arr, int l=0,int r=-1 ) {
	if ( r==-1 ) t_sz=r=arr.size()-1;
	Node *nd = create_node();
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
	Node *neo = create_node();
	if ( l >= r ) neo->dat= nd->dat+vl;
	else {
		int m =(l+r)/2;
		if ( x <= m ) neo->l=update(x,vl,nd->l,l,m),neo->r=nd->r;
		else neo->r=update(x,vl,nd->r,m+1,r),neo->l=nd->l;
		neo->dat = neo->l->dat + neo->r->dat;
	}
	return neo;
}

Node *bulk_update(const vector<pair<int,Data>>&vl,Node*nd ){
	// I assume elements of the form vl[i]={ind,x} where im supposed to
	// do Arr[vl[i].first]=vl[i].second (or whatever update im using
	int ind=0; //I assume vl is sorted by .first (index)
	function<Node*(Node*,int,int)> upd=[&](Node*nd, int l,int r){
		if ( ind==int(vl.size())|| r<vl[ind].first||vl[ind].first<l ) return nd;
		// debug(l,r,vl[ind].first,ind);
		Node *neo = create_node(); neo->dat=nd->dat;
		if ( l >= r ) while(ind<int(vl.size())&&vl[ind].first==l) neo->dat= neo->dat+vl[ind++].second;
		else {
			int m =(l+r)/2;
			neo->l=upd(nd->l,l,m); neo->r=upd(nd->r,m+1,r);
			neo->dat = neo->l->dat + neo->r->dat;
		}
		return neo;
	}; return upd(nd,0,t_sz);
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
/* EndOf persistent segtree */