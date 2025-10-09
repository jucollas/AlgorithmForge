/*
Autor: Oscar Vargas Pabon

No olvidar el guardar las nuevas versiones en rts[rt_n++]
*/

class Lazy{
public:
	int vl;
	Lazy(int v=0):vl(v){};
	
	void combine( const Lazy &l ){
		// unir con otro lazy (l estaba en el padre)
		vl+=l.vl+vl;
	}
	pair<Lazy,Lazy> separe( int l, int r ) {
		// separar en izquierda y derecha
		return {*this,*this};
	}
};
class Data {
	public:
	int sum;
	Data(int s=0) : sum(s) {};
	Data operator + ( const Data &o ) const {
		return Data(sum+o.sum);
	}
	void update( const Lazy &o, int l, int r ){
		// modificar el dato
		sum+=(r-l+1)*vl;
	}
};

class Node {
	public:
	Node *l,*r;
	Data dat;
	Lazy tag;
	Node(){ l=r=NULL; }
};

int t_sz, rt_n;
Node *rts[template_limit]; // roots;

list<Node*> to_erase; // guarda todo lo que se debe borrar;
Node* newNode(){ Node *neo = new Node(); to_erase.push_back(nd); return neo; }
void free_pseg(){
	for ( Node *nd : to_erase) delete nd;
	rt_n = 0; to_erase={};
}

Node* build( const vector<Data> &arr, int l=0,int r=-1 ) {
	if ( r==-1 ) t_sz=r=arr.size()-1;
	Node *nd = newNode();
	if ( l>=r ) nd->dat = arr[l];
	else {
		int m = (l+r)/2;
		nd->l = build(arr,l,m);
		nd->r = build(arr,m+1,r);
		nd->dat = nd->l->dat + nd->r->dat;
	}
	return nd;
}
void push( Node*nd, int l, int r ) {
/* Empuja el valor del tag a los nodos mas profundos */
	if ( nd->tag==Lazy() ) return;//no hay nada por pushear
	nd->dat.update( nd->tag,l,r );
	if ( l < r ) {
		Node *l=newNode(),*r=newNode();
		*l = *(nd->l); *r =*(nd->r);
		
		pair<Lazy,Lazy> tsep=nd->tag.separe(l,r);
		l->tag.combine(tsep.first);
		r->tag.combine(tsep.second);
		
		nd->l=l;nd->r=r;
	}
	nd->tag = Lazy(); // para no sobrecontar
}

Node *update( int ql,int qr, const Lazy &upd, Node *nd, int l=0, int r=t_sz ){
	push(nd,l,r);
	Node *neo = newNode();
	if ( ql<=l && r <= qr ) {
		neo=nd;
		neo->tag=upd;//funciona porque ya hize el update
		push(neo,l,r);
	} else if (!(r<ql||qr<l)){
		int m =(l+r)/2;
		neo->l=update(x,vl,nd->l,l,m);
		neo->r=update(x,vl,nd->r,m+1,r);
		
		neo->dat = neo->l->dat + neo->r->dat;
	} else neo = nd;
	return neo;
}
Data query( int ql, int qr, Node *nd, int l=0, int r=t_sz ) {
	Data res; push(nd,l,r);
	if ( ql<=l && r <= qr ) res = nd->dat;
	else if ( !(r<ql||qr<l) ) {
		int m = (l+r)/2;
		res = query(ql,qr,nd->l,l,m)+query(ql,qr,nd->r,m+1,r);
	}
	return res;
}

